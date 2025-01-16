set.seed(2024)
nsize = c(100, 200, 500, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6)
data_list  = list()
for(i in 1:length(nsize)){
  data_list[[i]] = tgcc:::make.mixgaussian(nsize[i])$data
}
saveRDS(data_list, file = "./simulation/Figure_4/datasets.rds")
