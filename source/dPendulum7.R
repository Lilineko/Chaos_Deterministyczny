library(tidyverse)
library(gganimate)
library(ggraph)

# constants
gravitationalConstant <- 9.81  # acceleration due to gravity [m/s^2]
p1Length <- 1.0                # length of pendulum 1 [m]
p1Mass <- 1.0                  # mass of pendulum 1 [kg]
p2Length <- 1.0                # length of pendulum 2 [m]
p2Mass <- 1.0                  # mass of pendulum 2 [kg]

# initial conditions
p1Angle <- rep(120.0, 7)               
p1AngularVelocity <- rep(0.0, 7)     
p2Angle <- rep(1500 : 1506) / 10
p2AngularVelocity <- rep(0.0, 7)  
dt = 0.001                      
totalTime = 15.0              

p1Angle = p1Angle * pi / 180
p1AngularVelocity = p1AngularVelocity * pi / 180
p2Angle = p2Angle * pi / 180
p2AngularVelocity = p2AngularVelocity * pi / 180

# definition of the function calculating data points for the double pendulum problem
get7Pendulum <- function(p1Angle, p1AngularVelocity, p2Angle, p2AngularVelocity, dt, totalTime) {
  # declarations - BEGIN
  
  time = 0.0
  p1AngularAcceleration = 0.0
  p2AngularAcceleration = 0.0
  
  nSteps = floor(totalTime / dt)
  pendulum = matrix(nrow = nSteps, ncol = 29)
  # declarations - END
  
  pendulum[1, 1:7] = p1Angle
  pendulum[1, 8:14] = p1AngularVelocity
  pendulum[1, 15:21] = p2Angle
  pendulum[1, 22:28] = p2AngularVelocity
  pendulum[1, 29] = time
  
  for(it in 2:nSteps) {
    # update of time
    time = time + dt
    
    # calculation of angular acceleration (mathematical formula can be found at myphysicslab.com)
    p1Numerator1 = -gravitationalConstant * ((2.0 * p1Mass + p2Mass) * sin(p1Angle) + p2Mass * sin(p1Angle - 2.0 * p2Angle))
    p1Numerator2 = -2.0 * p2Mass * sin(p1Angle - p2Angle) * ((p2AngularVelocity ^ 2) * p2Length + (p1AngularVelocity ^ 2) * p1Length * cos(p1Angle - p2Angle))
    p1Denominator = p1Length * (2.0 * p1Mass + p2Mass * (1.0 - cos(2.0 * (p1Angle - p2Angle))))
    
    p1AngularAcceleration = (p1Numerator1 + p1Numerator2) / p1Denominator
    
    p2Numerator1 = (p1Mass + p2Mass) * (gravitationalConstant * cos(p1Angle) + p1Length * (p1AngularVelocity ^ 2))
    p2Numerator2 = ((p2AngularVelocity ^ 2) * p2Length * p2Mass * cos(p1Angle - p2Angle))
    p2Denominator = p1Denominator * p2Length / p1Length
    
    p2AngularAcceleration = 2.0 * sin(p1Angle - p2Angle) * (p2Numerator1 + p2Numerator2) / p2Denominator
    
    # calculation of angular velocity
    p1AngularVelocity = p1AngularVelocity + p1AngularAcceleration * dt
    p2AngularVelocity = p2AngularVelocity + p2AngularAcceleration * dt
    
    # calculation of angles
    p1Angle = p1Angle + p1AngularVelocity * dt
    p2Angle = p2Angle + p2AngularVelocity * dt
    
    # storage of the data
    pendulum[it, 1:7] = p1Angle
    pendulum[it, 8:14] = p1AngularVelocity
    pendulum[it, 15:21] = p2Angle
    pendulum[it, 22:28] = p2AngularVelocity
    pendulum[it, 29] = time
  }
  
  return(pendulum)
}

# calculation of data points for the double pendulum problem
pendulum <- get7Pendulum(p1Angle, p1AngularVelocity, p2Angle, p2AngularVelocity, dt, totalTime)

# definition of the function converting radial coordinates to cartesian ones
getCartesianCoordinates <- function(radialCoordinates) {
  p1Angles = radialCoordinates[ , 1:7]
  p2Angles = radialCoordinates[ , 8:14]
  
  p1x = p1Length * sin(p1Angles)
  p1y = -p1Length * cos(p1Angles)
  p2x = p1x + p2Length * sin(p2Angles)
  p2y = p1y - p2Length * cos(p2Angles)
  
  cartesianCoordinates <- matrix(nrow = nrow(p1Angles), ncol = 29)
  cartesianCoordinates[ , 1:7] = p1x
  cartesianCoordinates[ , 8:14] = p1y
  cartesianCoordinates[ , 15:21] = p2x
  cartesianCoordinates[ , 22:28] = p2y
  cartesianCoordinates[ , 29] = radialCoordinates[ , 15]
  
  return(cartesianCoordinates)
}

coordinates <- getCartesianCoordinates(pendulum[ ,c(1:7, 15:21, 29)])

# animation
reduction = 10
reducedLength = dim(coordinates)[1] / reduction
time <- rep(0.0, reducedLength)
p1x <- matrix(nrow = reducedLength, ncol = 7)
p1y <- matrix(nrow = reducedLength, ncol = 7)
p2x <- matrix(nrow = reducedLength, ncol = 7)
p2y <- matrix(nrow = reducedLength, ncol = 7)
accumulatedTime <- time
for(it in 1 : reducedLength) {
  which = 1 + (it - 1) * reduction
  time[it] <- coordinates[which, 29]
  p1x[it, ] <- coordinates[which, 1:7]
  p1y[it, ] <- coordinates[which, 8:14]
  p2x[it, ] <- coordinates[which, 15:21]
  p2y[it, ] <- coordinates[which, 22:28]
}

drawingData <- tibble(time, p1x, p1y, p2x, p2y, group = 1)

ggplot(drawingData) +
  geom_segment(aes(x = -6.4, xend = 6.4, y = -2.4, yend = -2.4)) +
  geom_segment(aes(x = -6.4, xend = -6.4, y = -2.4, yend = 2.4)) +
  geom_segment(aes(x = 6.4, xend = 6.4, y = -2.4, yend = 2.4)) +
  geom_segment(aes(x = -6.4, xend = 6.4, y = 2.4, yend = 2.4)) +
  geom_segment(aes(x = 0, xend = 0, y = -2.4, yend = 2.4)) +
  geom_segment(aes(x = -3.2, y = 0, xend = p1x[ , 1] - 3.2, yend = p1y[ , 1])) +
  geom_segment(aes(x = -3.2, y = 0, xend = p1x[ , 2] - 3.2, yend = p1y[ , 2])) +
  geom_segment(aes(x = -3.2, y = 0, xend = p1x[ , 3] - 3.2, yend = p1y[ , 3])) +
  geom_segment(aes(x = -3.2, y = 0, xend = p1x[ , 4] - 3.2, yend = p1y[ , 4])) +
  geom_segment(aes(x = -3.2, y = 0, xend = p1x[ , 5] - 3.2, yend = p1y[ , 5])) +
  geom_segment(aes(x = -3.2, y = 0, xend = p1x[ , 6] - 3.2, yend = p1y[ , 6])) +
  geom_segment(aes(x = -3.2, y = 0, xend = p1x[ , 7] - 3.2, yend = p1y[ , 7])) +
  geom_segment(aes(x = p1x[ , 1] - 3.2, y = p1y[ , 1], xend = p2x[ , 1] - 3.2, yend = p2y[ , 1])) +
  geom_segment(aes(x = p1x[ , 2] - 3.2, y = p1y[ , 2], xend = p2x[ , 2] - 3.2, yend = p2y[ , 2])) +
  geom_segment(aes(x = p1x[ , 3] - 3.2, y = p1y[ , 3], xend = p2x[ , 3] - 3.2, yend = p2y[ , 3])) +
  geom_segment(aes(x = p1x[ , 4] - 3.2, y = p1y[ , 4], xend = p2x[ , 4] - 3.2, yend = p2y[ , 4])) +
  geom_segment(aes(x = p1x[ , 5] - 3.2, y = p1y[ , 5], xend = p2x[ , 5] - 3.2, yend = p2y[ , 5])) +
  geom_segment(aes(x = p1x[ , 6] - 3.2, y = p1y[ , 6], xend = p2x[ , 6] - 3.2, yend = p2y[ , 6])) +
  geom_segment(aes(x = p1x[ , 7] - 3.2, y = p1y[ , 7], xend = p2x[ , 7] - 3.2, yend = p2y[ , 7])) +
  geom_point(x = -3.2, y = 0, size  = 2.0) +
  geom_point(aes(p1x[ , 1] - 3.2, p1y[ , 1]), col = "dark gray", size = 5.0 * p1Mass) +
  geom_point(aes(p1x[ , 2] - 3.2, p1y[ , 2]), col = "dark gray", size = 5.0 * p1Mass) +
  geom_point(aes(p1x[ , 3] - 3.2, p1y[ , 3]), col = "dark gray", size = 5.0 * p1Mass) +
  geom_point(aes(p1x[ , 4] - 3.2, p1y[ , 4]), col = "dark gray", size = 5.0 * p1Mass) +
  geom_point(aes(p1x[ , 5] - 3.2, p1y[ , 5]), col = "dark gray", size = 5.0 * p1Mass) +
  geom_point(aes(p1x[ , 6] - 3.2, p1y[ , 6]), col = "dark gray", size = 5.0 * p1Mass) +
  geom_point(aes(p1x[ , 7] - 3.2, p1y[ , 7]), col = "dark gray", size = 5.0 * p1Mass) +
  geom_point(aes(p2x[ , 1] - 3.2, p2y[ , 1]), col = "magenta", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 2] - 3.2, p2y[ , 2]), col = "purple", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 3] - 3.2, p2y[ , 3]), col = "blue", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 4] - 3.2, p2y[ , 4]), col = "green", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 5] - 3.2, p2y[ , 5]), col = "yellow", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 6] - 3.2, p2y[ , 6]), col = "orange", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 7] - 3.2, p2y[ , 7]), col = "red", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 1] + 3.2, p2y[ , 1]), col = "magenta", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 2] + 3.2, p2y[ , 2]), col = "purple", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 3] + 3.2, p2y[ , 3]), col = "blue", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 4] + 3.2, p2y[ , 4]), col = "green", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 5] + 3.2, p2y[ , 5]), col = "yellow", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 6] + 3.2, p2y[ , 6]), col = "orange", size = 5.0 * p2Mass) +
  geom_point(aes(p2x[ , 7] + 3.2, p2y[ , 7]), col = "red", size = 5.0 * p2Mass) +
  geom_path(aes(p2x[ , 1] + 3.2, p2y[ , 1]), col = "magenta", alpha = 0.2) +
  geom_path(aes(p2x[ , 2] + 3.2, p2y[ , 2]), col = "purple", alpha = 0.2) +
  geom_path(aes(p2x[ , 3] + 3.2, p2y[ , 3]), col = "blue", alpha = 0.2) +
  geom_path(aes(p2x[ , 4] + 3.2, p2y[ , 4]), col = "green", alpha = 0.2) +
  geom_path(aes(p2x[ , 5] + 3.2, p2y[ , 5]), col = "yellow", alpha = 0.2) +
  geom_path(aes(p2x[ , 6] + 3.2, p2y[ , 6]), col = "orange", alpha = 0.2) +
  geom_path(aes(p2x[ , 7] + 3.2, p2y[ , 7]), col = "red", alpha = 0.2) +
  scale_y_continuous(limits = c(-2.4, 2.4)) +
  scale_x_continuous(limits = c(-6.4, 6.4)) +
  ggraph::theme_graph() +
  labs(title="{floor(frame_along)} s") +
  transition_reveal(time) -> pendulum7Plots
  pendulum7Plots <- pendulum7Plots + coord_fixed()

# pendulum7Animation <- animate(pendulum7Plots, nframes = nrow(drawingData) / 3, fps = 25)
# pendulum7Animation

animate(pendulum7Plots, nframes = nrow(drawingData) / 3, fps = 33)
anim_save("dPendulum7.mp4")

system("ffmpeg -i dPendulum7.mp4 -vf 'setpts=1*PTS' dPendulum7_f.mp4")
