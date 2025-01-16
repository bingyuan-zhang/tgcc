library(tgcc)
library(ggplot2)
set.seed(2024)
dl <- tgcc:::make.multidim.threeclust(n = 500, p = 3)

library(dplyr)
library(plotly)
colors <- c("red", "green", "blue")
point_colors <- colors[dl$label]

fig <-
  plotly::plot_ly(
    x = ~ dl$data[, 1],
    y = ~ dl$data[, 2],
    z = ~ dl$data[, 3],
    type = 'scatter3d',
    mode = 'markers',
    marker = list(color = point_colors, size = 1.5)
  )
fig %>% layout(scene = list(
  xaxis = list(title = ""),
  yaxis = list(title = ""),
  zaxis = list(title = "")
))

