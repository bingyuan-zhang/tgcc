# show result table
AC.mat <- biRes$AC;
ARI.mat <- biRes$ARI;
Time.mat <- biRes$Time;
################################################################################
# make the matrix
rowMeans(AC.mat[1,,]) %>% round(.,digits=3)
apply(AC.mat[1,,], 1, sd) %>% round(.,digits=3)

rowMeans(ARI.mat[1,,]) %>% round(.,digits=3)
apply(ARI.mat[1,,], 1, sd) %>% round(.,digits=3)

rowMeans(Time.mat[1,,]) %>% round(.,digits=3)
apply(Time.mat[1,,], 1, sd) %>% round(.,digits=3)


rowMeans(AC.mat[2,,]) %>% round(.,digits=3)
apply(AC.mat[2,,], 1, sd) %>% round(.,digits=3)

rowMeans(ARI.mat[2,,]) %>% round(.,digits=3)
apply(ARI.mat[2,,], 1, sd) %>% round(.,digits=3)

rowMeans(Time.mat[2,,]) %>% round(.,digits=3)
apply(Time.mat[2,,], 1, sd) %>% round(.,digits=3)
################################################################################
