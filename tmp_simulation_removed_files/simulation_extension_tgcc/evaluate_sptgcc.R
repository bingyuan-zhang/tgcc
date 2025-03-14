# sptgcc evaluation
source("./simulation/simulation_extension_tgcc/utils_sptgcc.R")
set.seed(2024)
gdata_list <- generate_spcc_datalist()
spRes <- evaluate_sptgcc(gdata_list)

# saveRDS(spRes, file="./spRes.rds")
