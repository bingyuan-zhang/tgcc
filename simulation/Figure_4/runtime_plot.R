# plot out
rt_dbscan <- readRDS("./simulation/Figure_4/rt_dbscan.rds")
rt_sc <- readRDS("./simulation/Figure_4/rt_sc.rds")
rt_km <- readRDS("./simulation/Figure_4/rt_km.rds")
rt_slc <- readRDS("./simulation/Figure_4/rt_slc.rds")
rt_carp <- readRDS("./simulation/Figure_4/rt_carp.rds")
rt_cpaint <- readRDS("./simulation/Figure_4/rt_cpaint.rds")
rt_ccmm <- readRDS("./simulation/Figure_4/rt_ccmm.rds")
rt_knn <- readRDS("./simulation/Figure_4/rt_knn.rds")
rt_mst <- readRDS("./simulation/Figure_4/rt_mst.rds")
rt_tgcc <- readRDS("./simulation/Figure_4/rt_tgcc.rds")

df <- rbind(rt_dbscan, rt_sc, rt_km, rt_slc,
  rt_carp, rt_cpaint, rt_ccmm,
  rt_knn, rt_mst, rt_tgcc)

df$rtime <- log10(df$rtime)
df$sample <- log10(df$sample)


methods <- c(
  "MST", "KNN",
  "DBS", "SC", "KM", "SLC",
  "CARP", "CPAINT", "CCMM", "TGCC"
)
col_methods <- c(
  "MST" = "#000000",
  "KNN" = "#999999",
  "DBS" = "#4DAF4A",
  "SC" = "#984EA3",
  "KM" = "#FF7F00",
  "SLC" = "#0000CD",
  "CARP" = "#A65628",
  "CPAINT" = "#F781BF",
  "CCMM" = "#377EB8",
  "TGCC" = "#E41A1C"
)
lty_methods <- c(
  "MST" = 3, "KNN" = 3,
  "DBS" = 2, "SC" = 2, "KM" = 2, "SLC" = 2,
  "CARP" = 1, "CPAINT" = 1, "CCMM" = 1, "TGCC" = 1
)
xlabels <- parse(text = paste0("10^", 2:6))
ylabels <- parse(text = paste0("10^", -4:2))
plot(NULL, type = "l", xlim = c(2.5, 6), ylim = c(-4, 2),
  xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(1, at = 2:6, labels = xlabels, cex.axis = 1.5)
axis(2, at = -4:2, labels = ylabels, cex.axis = 1.5, las = 2)

for(i in seq_along(methods)) {
  method <- methods[i]
  lty <- lty_methods[i]
  col <- col_methods[i]
  cur_df <- df[df$method == method,]
  lines(cur_df$sample, cur_df$rtime, lty = lty, col = col, lwd=2)
}
legend("bottomright", legend = methods, lty = lty_methods, col = col_methods, bty = "n", lwd = 2, ncol=2, cex = 1.2)

