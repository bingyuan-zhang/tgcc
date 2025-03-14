# sptgcc evaluation
source("./simulation/Table2-3/utils_sptgcc.R")
set.seed(2025)
gdata_list <- generate_spcc_datalist()
spRes <- evaluate_sptgcc(gdata_list)
saveRDS(spRes, file = "./simulation/Table2-3/spRes.rds")
