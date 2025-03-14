library(tgcc)
library(dendextend)

n <- 400
set.seed(2024)
dataset <- tgcc:::make.twomoons(n)
data <- dataset$data
color <- rep("red",n)
color[which(dataset$label == 1)] <- "blue"
plot(data, pch = 16, asp = 1, col = color)

# SLC
slc.twomoons <- hclust(d = dist(data), method = "single")
dend.slc.twomoons <- as.dendrogram(slc.twomoons)
labels_colors(dend.slc.twomoons) <- color[slc.twomoons$order]
dend <- dend.slc.twomoons
plot(dend.slc.twomoons)
colored_bars(color, dend)

# TGCC
res_lam <- tgcc:::checkLamTGCC(data, bandwidth = 10)
lambdaSeq <- seq(0.1, res_lam$lam_max, length.out = 100)
tgccFit <- tgcc::tgCC(data, lambdaSeq, bandwidth = 10)
pars <- tgcc:::parDendrogram(tgccFit)

h.tgcc.twomoons <- list(merge = pars$m, height = pars$h, order = pars$iorder)
class(h.tgcc.twomoons) <- "hclust"
dend.tgcc.twomoons <- as.dendrogram(h.tgcc.twomoons)
labels_colors(dend.tgcc.twomoons) <- color[h.tgcc.twomoons$order]
plot(dend.tgcc.twomoons)
colored_bars(color, dend.tgcc.twomoons)

# Figure 1
lab = tgcc::clusterLabel(tgccFit, numClusters = 2)
plot(data, col = c("red", "blue")[lab], asp = 1, pch = 16, xaxt = "n", yaxt = "n", ann=FALSE, cex = 0.7)
plot(dend.slc.twomoons, main = "SCL")
colored_bars(color, dend, add = TRUE)
plot(dend.tgcc.twomoons, main = "TGCC")
colored_bars(color, dend.tgcc.twomoons, add = TRUE)
