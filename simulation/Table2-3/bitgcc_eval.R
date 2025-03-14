# bitgcc evaluation
source("./simulation/Table2-3/utils_bitgcc.R")
set.seed(2025)
gdata_list <- generate_bicc_datalist()
biRes <- evaluate_bitgcc(gdata_list)
saveRDS(biRes, file = "./simulation/Table2-3/biRes.rds")

