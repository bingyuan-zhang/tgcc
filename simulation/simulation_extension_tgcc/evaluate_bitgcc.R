# bitgcc evaluation
source("./simulation/simulation_extension_tgcc/utils_bitgcc.R")
set.seed(2024)
gdata_list <- generate_bicc_datalist()
biRes <- evaluate_bitgcc(gdata_list)

# saveRDS(biRes, file="./biRes.rds")

