library(limma)
library(BayesMP)
## read in the data

data_brown <- read.csv("https://bayesmp.github.io/data/mouseMetabolism/data_brown.csv", row.names = 1)
data_heart <- read.csv("https://bayesmp.github.io/data/mouseMetabolism/data_heart.csv", row.names = 1)
data_liver <- read.csv("https://bayesmp.github.io/data/mouseMetabolism/data_liver.csv", row.names = 1)

## Verify gene names match across three tissues
all(rownames(data_brown) == rownames(data_heart))
all(rownames(data_brown) == rownames(data_liver))

## Combime these three studies as list
dataExp <- list(brown=data_brown, heart=data_heart, liver=data_liver)

## Check the dimension of the three studies
sapply(dataExp, dim)

## Check the head of the three studies
sapply(dataExp, head)

## perform differential expression analysis for each of these three tissues.

## Create an empty matrix to store Z value. Each row represents a gene and each column represent a study/tissue. Note that Z value matrix is the input for BayesMP method. Z value can be calculated by inverse CDF of the p-value from differential expression analysis. A positive Z value indicates that the gene is upregulated in that study/tissue.

Z <- matrix(0,nrow=nrow(dataExp[[1]]),ncol=length(dataExp))
rownames(Z) <- rownames(dataExp[[1]])
colnames(Z) <- names(dataExp)

for(s in 1:length(dataExp)){
	adata <- dataExp[[s]]
	ControlLabel = grep('wt',colnames(adata))
	caseLabel = grep('LCAD',colnames(adata))
	label <- rep(NA, ncol(adata))
	label[ControlLabel] = 0
	label[caseLabel] = 1
	
	design = model.matrix(~label)	## design matrix
	fit <- lmFit(adata,design)		## fit limma model
	fit <- eBayes(fit)		
	
	aeffectsize <- fit$coefficients[,2]	## get effect sizes
	Z[aeffectsize>0,s] <- -qnorm(fit$p.value[aeffectsize>0,2]/2)
	Z[aeffectsize<0,s] <- qnorm(fit$p.value[aeffectsize<0,2]/2)
	
}

## checkout top rows of Z
head(Z)


## Perform MCMC of BayesMP. There are two models the alternative distributions: Dirichlet Process (DP) or Gaussian distribution (DP). For the purpose of fast computing, we adopt DP method.
WD <- '~/Desktop/'
setwd(WD) ## You can set the working directory here. Some MCMC results will be saved here.

niter <- 1000 ## number of iterations.
burnin <- 200 ## number of burnin samples. Note the burnin samples will be disgarded.

set.seed(15213)	## set random seed

## writeHSall=T will save the intermediate results for the Bayesian hypothesis testing settings.
## writeY=T will save the intermediate results for the posterior samples of Y, which is used as the input for metaPattern detection.
## If you want to track the MCMC samples for other paramters including gamma, Pi, Delta, 
## One can save these results by set writeGamma, writePi, writeDelta to be TRUE.
## One can set gamma to be fixed by setting updateGamma = TRUE, if he has a good estimate of initial gamma.

system.time(BayesMP_DP(Z,niter=niter, burnin=burnin, writeY = T, writeHSall=T))



## Task I, meta analysis. There are three hypothesis testing settings.
HSallRes <- read.table('BayesMP_DP_HSall.txt')
## 
## HSB
HSb_belief <- HSallRes[,1]/(niter - burnin)
HSb_qvalue <- BayesianFDR(HSb_belief)
sum(HSb_qvalue<0.05)
sum(HSb_qvalue<0.01)


## HSA
S <- ncol(Z)
HSa_belief <- HSallRes[,S]/(niter - burnin)
HSa_qvalue <- BayesianFDR(HSa_belief)
sum(HSa_qvalue<0.05)

## HSr
r <- 2
HSr_belief <- HSallRes[,r]/(niter - burnin)
HSr_qvalue <- BayesianFDR(HSr_belief)
sum(HSr_qvalue<0.05)
 


## Task II, differential expression meta-pattern. 
## MetaPattern
con  <- file('BayesMP_DP_Y.txt', open = "r")

G <- nrow(Z)
S <- ncol(Z)

## the resYplus and resYminus matrices are used to save "latex pression"
resYplus <- matrix(0,G,S)
resYminus <- matrix(0,G,S)


i = 1
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  if(i>burnin){
	  #print(i)
	  seven = strsplit(oneLine, "\t")[[1]]
	  thisY <- matrix(as.numeric(seven),G,S)
  	
	  ## for individual studies
	  resYplus[thisY>0] <- resYplus[thisY>0] + 1
	  resYminus[thisY<0] <- resYminus[thisY<0] + 1
  }    
  i = i + 1
} 

close(con)
	
## we only consider HSb_qvalue<0.01 for the meta-pattern detection
resYplus_DE <- resYplus[HSb_qvalue<0.01,]
resYminus_DE <- resYminus[HSb_qvalue<0.01,]

## tight clustering
dissimilarity <- distance(resYplus_DE, resYminus_DE, niter - burnin) ## calculate dissimilarity matrix for each pair of genes.
atarget <- 2
## here, may need to depend on the consensus clustering algorithm to detect gene modules.
tightClustResult <- tightClustPam(dissimilarity, target=atarget, k.min=25) ## perform the tight clustering to identify meta-pattern modules. ## k.min controls the size of the final meta module. A small k.min will result in large module size and a large k.min will result in small module size.

## visualize the posterior probablity

par(mfrow = c(2,2))
for(i in 1:atarget){
	clustSelection <- tightClustResult$cluster == i
	numS = ncol(resYplus_DE)
	barData <- cbind(resYplus_DE,resYminus_DE)[clustSelection,]/(niter - burnin)
	barMean <- colMeans(barData)
	barSd <- apply(barData,2,sd)

	mp <- barplot(barMean, axes=FALSE, axisnames=FALSE, ylim=c(0, 1.2),
	              col=c(rep('red',numS), rep('blue',numS)), main=paste('target',i,'\n n =',sum(clustSelection)), 
				  xlab="Study", ylab="prosterior probablity")
	axis(1, labels=c(paste(1:numS,"+",sep=''), paste(1:numS,"-",sep='')), at = mp)
	axis(2, at=seq(0 , 1, by=0.2))
	 
	box()

	# Plot the vertical lines of the error bars
	# The vertical bars are plotted at the midpoints
	segments(mp, barMean - barSd, mp, barMean + barSd, lwd=2)
	# Now plot the horizontal bounds for the error bars
	# 1. The lower bar
	segments(mp - 0.1, barMean - barSd, mp + 0.1, barMean - barSd, lwd=2)
	# 2. The upper bar
	segments(mp - 0.1, barMean + barSd, mp + 0.1, barMean + barSd, lwd=2)
}

## visualize the meta-pattern for the first detected module

library(gplots)

for(s in 1:length(dataExp)){
	adata <- dataExp[[s]]
	aname <- names(dataExp)[s]
	bdata <- adata[HSb_qvalue<0.01, ][tightClustResult$cluster == 1 ,]
	cdata <- as.matrix(bdata)
	ddata <- t(scale(t(cdata))) ## standardize the data such that for each gene, the mean is 0 and sd is 1.

	ColSideColors <- rep("black", ncol(adata))
	ColSideColors[grep('LCAD',colnames(adata))] <- "red"
	
	heatmap.2(ddata,Rowv = NA,Colv = NA,
	          na.color="grey", labRow=FALSE, main = aname, 
	          col = "greenred", trace = "none", ColSideColors = ColSideColors)
}



