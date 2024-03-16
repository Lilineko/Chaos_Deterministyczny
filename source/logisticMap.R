library(tidyverse)
library(gganimate)
library(ggraph)

logisticMap <- function(lambda, n = 1000, m = 500) {
  xold = 0.01
  for (it in 1:n) {
    xnew = lambda * xold * (1.0 - xold)
    xold = xnew
  }
  result = matrix(nrow = m, ncol = 2)
  for(it in 1:m) {
    xnew = lambda * xold * (1.0 - xold)
    xold = xnew
    result[it, 1] = lambda
    result[it, 2] = xold
  }
  return(result)
}

min = 2.4
max = 4.0
kmax = 1600
m = 400
lambda = rep(0.0, (kmax + 1) * m)
attractor = rep(0.0, (kmax + 1) * m)
id = 1
for(ik in 0:kmax) {
  bif = logisticMap(min + (max - min) * (ik / kmax), 1000, m)
  lambda[id:(id + m - 1)] = bif[ , 1]
  attractor[id:(id + m - 1)] = bif[ , 2]
  id = id + m
}

drawingData <- tibble(lambda, attractor, group = 1)

ggplot(drawingData) +
  geom_segment(aes(x = min, xend = max, y = 0, yend = 0)) +
  geom_segment(aes(x = min, xend = max, y = 1, yend = 1)) +
  geom_segment(aes(x = min, xend = min, y = 0, yend = 1)) +
  geom_segment(aes(x = max, xend = max, y = 0, yend = 1)) +
  geom_point(aes(lambda, attractor), col = "red", size = 0.5, alpha = 0.05) +
  scale_y_continuous(limits = c(0.0, 1.0)) +
  scale_x_continuous(limits = c(min, max)) +
  labs(x = "r", y = "x") +
  theme(aspect.ratio = 0.618, panel.background = element_blank()) +
  shadow_mark(size = 0.5, alpha = 0.02, color = "blue", future = TRUE) +
  transition_time(lambda) -> bifurcationPlots

animate(bifurcationPlots, nframes = 100, fps = 25, width = 1000, height = 618)
anim_save(filename = "bifurcationDiagram.gif")
