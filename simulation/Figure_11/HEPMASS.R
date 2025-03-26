library(readr)
library(tgcc)
HEPMASS <- read_csv("../newrealdata/1000_test.csv.gz")
L = HEPMASS[,1]
L = L$`# label`

Alldata = HEPMASS[,-1]
data <- as.matrix(Alldata[1:1e6,])
rdata <- scale(data[1:1e6, c(7, 11, 15, 26, 27)])
rlabel <- L[1:1e6]

rm(L)
rm(Alldata)
rm(HEPMASS)

lamseq <- seq(from = 100, to = nrow(data)/2, length.out = 1000)
stime <- proc.time()
tgcc.fit <- tgCC(data = rdata, lambdaSeq = lamseq, bandwidth = 50, probThresh = 0.01)
etime <- proc.time()
etime - stime
estlabel <- clusterLabel(tgcc.fit, numClusters = 2)
ari <- mclust::adjustedRandIndex(estlabel, rlabel)
er <- mclust::classError(estlabel, rlabel)$errorRate
1 - er; ari

source("./simulation/Figure_11/utils.R")
res.heat <- parsHeatmap(tgcc.fit, rlabel, show = 500, k = 2)
datanew <- res.heat$datanew
dend <- res.heat$dend
labelnew <- res.heat$rowcolor
col_new <- rep("orangered",nrow(datanew))
col_new[which(labelnew==0)] <- "royalblue"
heatmap(
  datanew,
  Rowv = dend,
  Colv = NA,
  xlab = "",
  labRow = FALSE,
  labCol = FALSE,
  RowSideColors = col_new,
  col = gplots::redgreen(500)
)
