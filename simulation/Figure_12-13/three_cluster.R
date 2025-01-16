set.seed(2024)
n <- 400
dataset <- tgcc:::make.threeclust(n)
color <- rep("red",n)
color[which(dataset$label == 1)] <- "blue"
color[which(dataset$label == 2)] <- "green"

# Figure 12
plot(dataset$data, pch = 16, asp = 1, col = color, xlab = "", ylab = "")
# add outlier points
outliers1 <- c(-3, -4.5)
outliers2 <- c(6, -3)
points(rbind(outliers1, outliers2), pch = 24, bg = "gold", col = "black", cex = 1.2)

# Figure 13
library(dendextend)
data <- rbind(dataset$data, outliers1, outliers2)
label <- c(dataset$label, "outlier", "outlier")
newcolor <- c(color, rep("gold", 2))

lambdaSeq <- seq(0.1, 250, length.out = 100)
tgccFit <- tgcc::tgCC(data, lambdaSeq, bandwidth = 10, depthThresh = 0, probThresh = 0)

pars <- tgcc:::parDendrogram(tgccFit)
h.tgcc <- list(merge = pars$m, height = pars$h, order = pars$iorder)
class(h.tgcc) <- "hclust"
dend.tgcc <- as.dendrogram(h.tgcc)
labels_colors(dend.tgcc) <- newcolor[h.tgcc$order]
branch_color <- ifelse(newcolor[h.tgcc$order] == "gold", "gold", "black")
dend.tgcc <- color_branches(dend.tgcc, k=n+2, col=branch_color)
plot(dend.tgcc)
colored_bars(newcolor, dend.tgcc)

lambdaSeq <- seq(0.1, 250, length.out = 100)
tgccFit <- tgcc::tgCC(data, lambdaSeq, bandwidth = 10, depthThresh = 10, probThresh = 0.1)

pars <- tgcc:::parDendrogram(tgccFit)
h.tgcc <- list(merge = pars$m, height = pars$h, order = pars$iorder)
class(h.tgcc) <- "hclust"
dend.tgcc <- as.dendrogram(h.tgcc)
labels_colors(dend.tgcc) <- newcolor[h.tgcc$order]
branch_color <- ifelse(newcolor[h.tgcc$order] == "gold", "gold", "black")
dend.tgcc <- color_branches(dend.tgcc, k=n+2, col=branch_color)
plot(dend.tgcc)
colored_bars(newcolor, dend.tgcc)
