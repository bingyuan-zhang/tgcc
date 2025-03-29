## Tree-guided convex clustering

R package `tgcc` implements the tree-guided convex clustering (TGCC) and its extensions.

## Installation

Run the following code in Rstudio to install `tgcc` from GitHub.

``` r
library(devtools)
install_github("bingyuan-zhang/tgcc")
```

## TGCC Example

We demonstrate TGCC using the following examples.

At the beginning, we load the `tgcc` R package.

``` r
library(tgcc)
```

The first example is a mixture Gaussian model as follows.

The ground truth of labels were shown in different colors as follows.

``` r
set.seed(2024)
dl <- tgcc:::make.mixgaussian(n = 400)
data <- dl$data
label <- dl$label
plot(data, col=label)
```

<img src="./inst/example_png/mixtureGaussian.svg" width="80%" style="display: block; margin: auto;"/>

To apply TGCC, we set a tuning parameter sequence for lambda, and the bandwidth for Gaussian kernel.

`tgcc::makeDendrogram` provides a dendrogram visualization of the convex clustering result. The y-axis represents the value of lambda, and colors indicate the input labels.

``` r
lambdaSeq <- seq(1, 200, length.out = 100)
bandwidth <- 10
tgccFit <- tgCC(
  data = data, 
  lambdaSeq = lambdaSeq, 
  bandwidth = bandwidth)
tgcc::makeDendrogram(tgccFit, label)
```

<img src="./inst/example_png/dendrogram1.svg" width="80%" style="display: block; margin: auto;"/>

By comparing the obtained labels with the ground truth, we can obtain the clustering performance accuracy (see also `classError` in `mclust` package).

``` r
predLabel <- tgcc::clusterLabel(tgccFit, numClusters=3)
1-mclust::classError(label, predLabel)$errorRate
```

```         
[1] 0.9875
```

## Extensions of TGCC

TGCC can be easily extended for other different settings under the convex clustering framework.

We show the extensions of TGCC in the sparse clustering and the biclustering setting.

#### Sparse Clustering Setting

For the sparse clustering setting, we generate the Four Spherical example.

`spTGCC` obtains the clusters while removing the noisy features simultaneously.

``` r
# generate FS model ground truth
set.seed(2024)
FSmodel <- tgcc:::make.fourspherical(n=400)
order <- order(FSmodel$label)
data <- FSmodel$data[order, ]
label <- FSmodel$label[order]
data0 <- FSmodel$groundtruth[order, ]
```

The estimated centroid matrices is visualized as heatmaps.

``` r
# tuning parameters
gammaSeq <- c(1, 1.5, 1.9)
lambdaSeq <- c(30, 60, 120)

# fit the spTGCC model
tgccFit <- spTGCC(
  data = data, 
  lambdaSeq = lambdaSeq,
  gammaSeq = gammaSeq,
  threshold = 1e-05,
  maxIter = 100)

# set color
range <- quantile(data[!is.na(data)], c(0.0, 0.5, 1))
colfun <- circlize::colorRamp2(range, 
  c("blue", "white", "red"))
breaks <- seq(min(data, na.rm = TRUE), max(data, na.rm = TRUE), length.out = 256)
colvec <- colfun(breaks)

showHeatmap <- function(data, colvec) {
  pheatmap::pheatmap(
    mat = data,
    color = colvec,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = FALSE,
    na_col = "grey65"
  )
}

showHeatmap(data0, colvec)
showHeatmap(data, colvec)
showHeatmap(tgccFit$theta[[1]], colvec)
showHeatmap(tgccFit$theta[[2]], colvec)
showHeatmap(tgccFit$theta[[3]], colvec)
```

<img src="./inst/example_png/sptgcc.png" width="100%" style="display: block; margin: auto;"/>

#### Biclustering Setting

For the biclustering setting, we generate the Check Board example.

`biTGCC` obtains the clusters in both samples and features.

``` r
set.seed(2024)
CBmodel <- tgcc:::make.checkerboard(n = 400)
data <- CBmodel$data
label <- CBmodel$label
data0 <- CBmodel$groundTruth
```

The estimated centroid matrices is visualized as heatmaps.

``` r
# tuning parameter of lambda and gamma
lambdaSeq <- gammaSeq <- c(20, 50, 150)

# fit the biTGCC model
tgccFit <- biTGCC(
    data,
    lambdaSeq,
    gammaSeq,
    threshold = 1e-05 * nrow(data) * ncol(data),
    maxIter = 100)

# set color
range <- quantile(data[!is.na(data)], c(0, 0.5, 1))
colfun <- circlize::colorRamp2(range, c("blue", "white", "red"))
breaks <- seq(min(data, na.rm = TRUE), max(data, na.rm = TRUE), length.out = 256)
colvec <- colfun(breaks)

showHeatmap(data0, colvec)
showHeatmap(data, colvec)
showHeatmap(tgccFit$theta[[1]], colvec)
showHeatmap(tgccFit$theta[[2]], colvec)
showHeatmap(tgccFit$theta[[3]], colvec)
```

<img src="./inst/example_png/bitgcc.png" width="100%" style="display: block; margin: auto;"/>

## Reference

Tree-Guided $L_1$-Convex Clustering by Bingyuan Zhang, and Yoshikazu Terada.

See also the [simulation code](https://github.com/bingyuan-zhang/tgcc_simulation_code) in the TGCC paper.
