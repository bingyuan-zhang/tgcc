
label = partdata$`# label`
data = as.matrix(partdata[,2:28])
n <- 10^3
x = data[1:n,]
xlabel = label[1:n]
lamseq = 2^seq(0.1, 100, 1)
dim(x)
rm(data)
rm(partdata)
#HEPMASS_10000
#===========================================
lamseq = 2^seq(10,100, 10)#2^seq(0.1, 100, 1)
#f6, f10, f14, f25, f26
#x <- scale(x[,-c(6, 10, 14, 18, 22)])
x <- scale(x[,c(7, 11, 15, 26, 27)])
pairs(x, pch = 16, cex = 0.3, col = xlabel + 1)

