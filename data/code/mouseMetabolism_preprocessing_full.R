options(stringsAsFactors = F)

WD <- '~/Desktop/BayesMP.github.io/data/mouseMetabolism/'
setwd(WD)

load('~/Desktop/BayesMP.github.io/data/raw/mouse.Rdata')
annotationFile <- read.delim('~/Desktop/BayesMP.github.io/data/raw/GPL1261-56135.txt',skip=16)
annotationDF <- data.frame(probes=annotationFile$ID, genes=sapply(strsplit(annotationFile$Gene.Symbol,split=' /// '),function(x) x[1]) )




data_brown0 <- DList[[1]]
data_brown1 <- data_brown0[,grep("b.wt|b.LCAD", colnames(data_brown0))]
rowMeans_brown <- rowMeans(data_brown1)
rowRanks_brown <- rank(rowMeans_brown)

data_heart0 <- DList[[2]]
data_heart1 <- data_heart0[,grep("h.wt|h.LCAD", colnames(data_heart0))]
rowMeans_heart <- rowMeans(data_heart1)
rowRanks_heart <- rank(rowMeans_heart)

data_liver0 <- DList[[3]]
data_liver1 <- data_liver0[,grep("l.wt|l.LCAD", colnames(data_liver0))]
rowMeans_liver <- rowMeans(data_liver1)
rowRanks_liver <- rank(rowMeans_liver)

# sumRanks <- rowRanks_brown + rowRanks_heart + rowRanks_liver
# filterIndex <- sumRanks > median(sumRanks, 0.5)
# annotationDF_filter <- annotationDF[match(names(filterIndex)[filterIndex], annotationDF$probes), ]

# data_brown2 <- data_brown1[filterIndex, ]
# data_heart2 <- data_heart1[filterIndex, ]
# data_liver2 <- data_liver1[filterIndex, ]

# annotationDF_match <- annotationDF[match(names(filterIndex)[filterIndex], annotationDF$probes), ]

annotationDF_filter <- annotationDF[match(rownames(data_brown1), annotationDF$probes), ]

data_brown2 <- data_brown1
data_heart2 <- data_heart1
data_liver2 <- data_liver1

all(rownames(data_brown2) == annotationDF_filter$probes)
all(rownames(data_heart2) == annotationDF_filter$probes)
all(rownames(data_liver2) == annotationDF_filter$probes)

rownames(data_brown2) <- annotationDF_filter$genes
rownames(data_heart2) <- annotationDF_filter$genes
rownames(data_liver2) <- annotationDF_filter$genes

data_brown <- data_brown2[!duplicated(annotationDF_filter$genes) & !is.na(annotationDF_filter$genes),]
data_heart <- data_heart2[!duplicated(annotationDF_filter$genes) & !is.na(annotationDF_filter$genes),]
data_liver <- data_liver2[!duplicated(annotationDF_filter$genes) & !is.na(annotationDF_filter$genes),]

write.csv(data_brown, "data_brown.csv")
write.csv(data_heart, "data_heart.csv")
write.csv(data_liver, "data_liver.csv")

