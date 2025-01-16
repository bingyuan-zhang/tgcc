source("./simulation/runtime_by_dimension/multidimension_data_generator.R")
library(tgcc)
library(ggplot2)
# Generate data
n_vec <- c(10000,20000)
p_vec <- c(5,10,15,20,25,30)
Repeat <- 10
data_list <- list()

set.seed(123)
for(s in 1:length(n_vec)){
  n <- n_vec[s]
  for(t in 1:length(p_vec)){
    p <- p_vec[t]
    for(r in 1:Repeat){
      dl = generate_multidim(n = n, p = p)
      data_list[[paste(n)]][[paste(p)]][[paste(r)]] = dl
    }
  }
}

# save computational time
lambdaSeq <- seq(1, 100000, length.out = 100)
bandwidth <- 1
resDF <- data.frame()
for(n in 1:length(n_vec)) for(p in 1:length(p_vec)) {
  tgccRT <- 0
  mstRT <- 0
  for(r in 1:Repeat){
    print(c(n, p, r))

    dl <- data_list[[n]][[p]][[r]]
    data <- dl$data
    label <- dl$label

    t1 <- proc.time()
    tgccFit <- tgCC(
      data = data,
      lambdaSeq = lambdaSeq,
      bandwidth = bandwidth,
      isNaive = FALSE
    )
    t2 <- proc.time()

    tgccRT <- tgccRT + tgccFit$tgccTime
    mstRT <- mstRT + tgccFit$mstTime
  }
  newDF1 <-
    data.frame(
      n = n_vec[n],
      p = p_vec[p],
      t = tgccRT / Repeat,
      method = "tgcc fit"
    )
  newDF2 <-
    data.frame(
      n = n_vec[n],
      p = p_vec[p],
      t = mstRT / Repeat,
      method = "mst fit"
    )
  resDF <- rbind(resDF, newDF1, newDF2)
}

ggplot(resDF,
  aes(
    x = p,
    y = log(t),
    group = interaction(n, method),
    col = interaction(n, method)
  )) +
  geom_line() +
  labs(x = "p", y = "Run Time (second)", color = "n and Type") +
  theme_minimal() +
  theme(legend.position = "right")

set.seed(1234)
dl = generate_multidim(n = 500, p = 3)

# 绘制3D散点图
library(dplyr)
library(plotly)
colors <- c("red", "green", "blue") # 1 -> red, 2 -> green, 3 -> blue
point_colors <- colors[dl$label]
# 绘制3D散点图
fig <- plotly::plot_ly(x = ~dl$data[,1], y = ~dl$data[,2], z = ~dl$data[,3],
  type = 'scatter3d', mode = 'markers',
  marker = list(color = point_colors, size = 1.5))
fig %>% layout(scene = list(
  xaxis = list(title = ""),
  yaxis = list(title = ""),
  zaxis = list(title = "")
))

