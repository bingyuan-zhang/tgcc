library(tgcc)
library(dendextend)

n <- 400
set.seed(2024)

dataset <- tgcc:::make.mixgaussian(n)
color <- rep("red",n)
color[which(dataset$label == 1)] <- "blue"
color[which(dataset$label == 2)] <- "green"
plot(dataset$data, pch = 16, asp = 1, col = color, xlab = "", ylab = "")

dataset <- tgcc:::make.threeclust(n)
color <- rep("red",n)
color[which(dataset$label == 1)] <- "blue"
color[which(dataset$label == 2)] <- "green"
plot(dataset$data, pch = 16, asp = 1, col = color, xlab = "", ylab = "")

dataset <- tgcc:::make.twocircles(n)
color <- rep("red",n)
color[which(dataset$label == 1)] <- "blue"
plot(dataset$data, pch = 16, asp = 1, col = color, xlab = "", ylab = "")

