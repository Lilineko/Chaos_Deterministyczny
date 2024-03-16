library(tidyverse)
library(gganimate)
library(ggraph)

# initial conditions (we consider two systems with slightly different initial position)
x <- c(2.0, 2.01)
y <- c(4.0, 4.0)
z <- c(8.0, 8.0)

dt <- 0.005
totalTime <- 60.0

getLorenzAttractor <- function(x, y, z, dt, totalTime) {
  # declarations - BEGIN
  sigma = 10.0
  rho = 28.0
  beta = 8.0 / 3.0
  
  time = 0.0
  
  nSteps = floor(totalTime / dt)
  attractor = matrix(nrow = nSteps, ncol = 5)  
  # declarations - END
  
  attractor[1, 1] = time
  attractor[1, 2] = x[1]
  attractor[1, 3] = z[1]
  attractor[1, 4] = x[2]
  attractor[1, 5] = z[2]
  
  for(it in 2:nSteps) {
    # update of time
    time = time + dt
    
    # calculations of Lorenz system
    dx = sigma * c(y[1] - x[1], y[2] - x[2]) * dt
    dy = c(x[1] * (rho - z[1]) - y[1], x[2] * (rho - z[2]) - y[2]) * dt
    dz = c(x[1] * y[1] - beta * z[1], x[2] * y[2] - beta * z[2]) * dt
    
    x = x + dx
    y = y + dy
    z = z + dz
    
    attractor[it, 1] = time
    attractor[it, 2] = x[1]
    attractor[it, 3] = z[1]
    attractor[it, 4] = x[2]
    attractor[it, 5] = z[2]
  }
  
  return(attractor)
}

attractor <- getLorenzAttractor(x, y, z, dt, totalTime)

# animation
time <- attractor[1, 1]
a1x <- attractor[1, 2]
a1y <- attractor[1, 3]
a2x <- attractor[1, 4]
a2y <- attractor[1, 5]
accumulatedTime <- time
for(it in 2 : dim(attractor)[1]) {
    time[length(time) + 1] <- attractor[it, 1]
    a1x[length(a1x) + 1] <- attractor[it, 2]
    a1y[length(a1y) + 1] <- attractor[it, 3]
    a2x[length(a2x) + 1] <- attractor[it, 4]
    a2y[length(a2y) + 1] <- attractor[it, 5]
}

drawingData <- tibble(time, a1x, a1y, a2x, a2y, group = 1)

ggplot(drawingData) +
  geom_segment(aes(x = -30, xend = 30, y = 0, yend = 0)) +
  geom_segment(aes(x = -30, xend = -30, y = 0, yend = 60)) +
  geom_segment(aes(x = 30, xend = 30, y = 0, yend = 60)) +
  geom_segment(aes(x = -30, xend = 30, y = 60, yend = 60)) +
  geom_point(aes(a1x, a1y), col = "red", size = 3.0) +
  geom_path(x = a1x, y = a1y, color = "red", alpha = 0.2) + 
  geom_point(aes(a2x, a2y), col = "blue", size = 3.0) +
  geom_path(x = a2x, y = a2y, color = "blue", alpha = 0.2) + 
  scale_y_continuous(limits = c(0, 60)) +
  scale_x_continuous(limits = c(-30, 30)) +
  ggraph::theme_graph() +
  #labs(title="{floor(frame_along)} s") +
  #shadow_mark(size = 0.2, alpha = 0.5) +
  #transition_time(time) -> attractorPlots
  transition_reveal(time) -> attractorPlots
  attractorPlots
  
# attractorAnimation <- animate(attractorPlots, nframes = nrow(drawingData) / 5, fps = 25)
# attractorAnimation

# animate(attractorPlots, nframes = nrow(drawingData) / 5, fps = 25)
# anim_save("lorenz.mp4")
# 
# system("ffmpeg -i lorenz.mp4 -vf 'setpts=2*PTS' lorenz_slow.mp4")
