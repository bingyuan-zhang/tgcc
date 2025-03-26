library(readr)
library(mlbench)
library(rmatio)
library(openxlsx)
source("./simulation/Table4/utils_table4.R")

# wine, 178*13, K = 3
wine <- read_csv("../newrealdata/wine.data", col_names = FALSE)
wine <- as.matrix(wine)
wine_data <- as.matrix(wine[,-1])
wine_data <- scale(wine_data)
wine_label <- as.numeric(wine[,1])

res_wine <- run_methods(wine_data, wine_label)
saveRDS(res_wine, file = "./simulation/Table4/res_wine.rds")

# Breast Cancer, 449*9, K = 2
# install.packages("mlbench")
data("BreastCancer")
d <- BreastCancer[ , 2:11]
d <- d[!duplicated(d), ]
d <- d[complete.cases(d), ]
bc_label <- as.numeric(d$Class)
d <- as.matrix(d[ , 1:9])
bc_data <- apply(d, 2, as.numeric)
bc_data <- scale(bc_data)

res_bc <- run_methods(bc_data, bc_label)
saveRDS(res_bc, file = "./simulation/Table4/res_bc.rds")

# Segmentation, 2310*18, K = 7
segmentation1 <- read_csv("../newrealdata/segmentation.data", col_names = FALSE)
segmentation2 <- read_csv("../newrealdata/segmentation.test", col_names = FALSE)
seg <- as.matrix(rbind(segmentation1, segmentation2))
seg_label <- factor(seg[,1])
seg_data <- matrix(nrow = nrow(seg), ncol = ncol(seg)-2)
for(i in 1:nrow(seg_data)){
  seg_data[i,] = as.numeric(seg[i,-c(1,4)])
}
seg_data <- scale(seg_data)

methods <- c("slc", "clc", "sc", "km", "dbscan", "cpaint", "ccmm", "tgcc")
res_seg <- run_methods(seg_data, seg_label, methods)
saveRDS(res_seg, file = "./simulation/Table4/res_seg.rds")

# mnist (10000*10, K = 10)
partmnist <- read.mat("../newrealdata/mnist_10000.mat")
L <- t(partmnist$labels)
data <- t(partmnist$X)
mnist_label <- rep(0, 10000)
for(i in 1:nrow(L)) mnist_label[i] <- which(L[i,] != 0)
mnist_data <- scale(data)

methods <- c("slc", "clc", "sc", "km", "dbscan", "cpaint", "ccmm", "tgcc")
res_mnist <- run_methods(mnist_data, mnist_label, methods)
saveRDS(res_mnist, file = "./simulation/Table4/res_mnist.rds")

# pendigit (10992*16, K = 10)
pendigits <- read_csv("../newrealdata/pendigits_csv.csv")
pend_label <- as.vector(pendigits$class)
pend_data <- as.matrix(pendigits[,-17])
pend_data <- scale(pend_data)

methods <- c("slc", "clc", "sc", "km", "dbscan", "cpaint", "ccmm", "tgcc")
res_pend <- run_methods(pend_data, pend_label, methods)
saveRDS(res_pend, file = "./simulation/Table4/res_pend.rds")


# Drybean 13611 * 16, K = 7
# install.packages("openxlsx")
Dry_Bean_Dataset <- read.xlsx("../newrealdata/Dry_Bean_Dataset.xlsx", sheet = 1)
label <- Dry_Bean_Dataset$Class
unique_label <- unique(label)
cov1 = list()
for(i in 1:7) cov1[[unique_label[i]]] = i
newlabel = label
for(i in 1:length(label)){
  newlabel[i] = cov1[[label[i]]]
}
dryb_label <- as.numeric(newlabel)
dryb_data <- as.matrix(Dry_Bean_Dataset[,-17])
dryb_data <- scale(dryb_data)

methods <- c("slc", "clc", "sc", "km", "dbscan", "cpaint", "ccmm", "tgcc")
res_dryb <- run_methods(dryb_data, dryb_label, methods)
saveRDS(res_dryb, file = "./simulation/Table4/res_dryb.rds")


# HEPMASS small, 100000*5, K = 2
HEPMASS <- read_csv("../newrealdata/1000_test.csv.gz")
L = HEPMASS[,1]
L = L$`# label`
Alldata = HEPMASS[,-1]
data <- as.matrix(Alldata[1:100000,])
hep_small_data <- scale(data[, c(7, 11, 15, 26, 27)])
hep_small_label <- L[1:100000]

methods <- c("km", "tgcc")
res_hep_small_1 <- run_methods(hep_small_data, hep_small_label, methods)
methods <- c("ccmm")
res_hep_small_2 <- run_methods(hep_small_data, hep_small_label, methods)
saveRDS(res_hep_small_1, file = "./simulation/Table4/res_hep_small_1.rds")
saveRDS(res_hep_small_2, file = "./simulation/Table4/res_hep_small_2.rds")

# CovType
covtype <- read_csv("../newrealdata/covtype.data", col_names = FALSE)
covt_label <- covtype$X55
covt_data <- covtype[,1:11]
covt_data <- scale(covt_data)

methods <- c("km")
res_covt_1 <- run_methods(covt_data, covt_label, methods)
saveRDS(res_covt_1, file = "./simulation/Table4/res_covt_1.rds")

methods <- c("tgcc")
res_covt_2 <- run_methods(covt_data, covt_label, methods)
saveRDS(res_covt_2, file = "./simulation/Table4/res_covt_2.rds")



