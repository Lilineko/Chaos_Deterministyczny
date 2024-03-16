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
p1Angle <- 120.0               # angle of pendulum 1 [degree]
p1AngularVelocity <- 0.0       # angular velocity of pendulum 1
p2Angle <- 150.0                # angle of pendulum 2 [degree]
p2AngularVelocity <- 0.0       # angular velocity of pendulum 2
dt = 0.001                      # time step [s]
totalTime = 60.00               # total time of animation [s]

init <- c(p1Angle, p1AngularVelocity, p2Angle, p2AngularVelocity) * pi / 180 # [+ conversion to radians]

# definition of the function calculating data points for the double pendulum problem
getDoublePendulum <- function(init, dt, totalTime) {
  # declarations - BEGIN
  
  time = 0.0
  p1Angle = init[1]
  p1AngularVelocity = init[2]
  p1AngularAcceleration = 0.0
  p2Angle = init[3]
  p2AngularVelocity = init[4]
  p2AngularAcceleration = 0.0
  
  nSteps = floor(totalTime / dt)
  pendulum = matrix(nrow = nSteps, ncol = length(init) + 1)
  # declarations - END
  
  pendulum[1, 1] = time
  pendulum[1, 2] = p1Angle
  pendulum[1, 3] = p1AngularVelocity
  pendulum[1, 4] = p2Angle
  pendulum[1, 5] = p2AngularVelocity
    
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
    pendulum[it, 1] = time
    pendulum[it, 2] = p1Angle
    pendulum[it, 3] = p1AngularVelocity
    pendulum[it, 4] = p2Angle
    pendulum[it, 5] = p2AngularVelocity
    # note that not every step of calculations must became a frame of animation but for now we dont care
    # data will be reduced only in the function creating animation
  }
    
  return(pendulum)
}

# calculation of data points for the double pendulum problem
pendulum <- getDoublePendulum(init, dt, totalTime)

# definition of the function converting radial coordinates to cartesian ones
getCartesianCoordinates <- function(radialCoordinates) {
  p1Angles = radialCoordinates[ , 2]
  p2Angles = radialCoordinates[ , 3]
  
  p1x = p1Length * sin(p1Angles)
  p1y = -p1Length * cos(p1Angles)
  p2x = p1x + p2Length * sin(p2Angles)
  p2y = p1y - p2Length * cos(p2Angles)
  
  cartesianCoordinates <- matrix(nrow = length(p1Angles), ncol = 5)
  cartesianCoordinates[ , 1] = radialCoordinates[ , 1]
  cartesianCoordinates[ , 2] = p1x
  cartesianCoordinates[ , 3] = p1y
  cartesianCoordinates[ , 4] = p2x
  cartesianCoordinates[ , 5] = p2y
  
  return(cartesianCoordinates)
}

coordinates <- getCartesianCoordinates(pendulum[ ,c(1, 2, 4)])

# animation
reduction = 20
time <- coordinates[1, 1]
p1x <- coordinates[1, 2]
p1y <- coordinates[1, 3]
p2x <- coordinates[1, 4]
p2y <- coordinates[1, 5]
accumulatedTime <- time
for(it in 1 : (dim(coordinates)[1] / reduction - 1)) {
    which = 1 + it * reduction
    time[length(time) + 1] <- coordinates[which, 1]
    p1x[length(p1x) + 1] <- coordinates[which, 2]
    p1y[length(p1y) + 1] <- coordinates[which, 3]
    p2x[length(p2x) + 1] <- coordinates[which, 4]
    p2y[length(p2y) + 1] <- coordinates[which, 5]
}

drawingData <- tibble(time, p1x, p1y, p2x, p2y, group = 1)

ggplot(drawingData) +
  geom_segment(aes(x = -2.5, xend = 2.5, y = -2.5, yend = -2.5)) +
  geom_segment(aes(x = -2.5, xend = -2.5, y = -2.5, yend = 2.5)) +
  geom_segment(aes(x = 2.5, xend = 2.5, y = -2.5, yend = 2.5)) +
  geom_segment(aes(x = -2.5, xend = 2.5, y = 2.5, yend = 2.5)) +
  geom_segment(aes(x = 0, y = 0, xend = p1x, yend = p1y)) +
  geom_segment(aes(x = p1x, y = p1y, xend = p2x, yend = p2y)) +
  geom_point(x = 0, y = 0, size  = 2.0) +
  geom_point(aes(p1x, p1y), col = "red", size = 5.0 * p1Mass) +
  geom_point(aes(p2x, p2y), col = "blue", size = 5.0 * p2Mass) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  ggraph::theme_graph() +
  labs(title="{floor(2 * frame_time)} s") +
  shadow_wake(wake_length = 1.0 / totalTime, alpha = FALSE, wrap = FALSE) +
  transition_time(time) -> pendulumPlots

# pendulumAnimation <- animate(pendulumPlots, nframes = nrow(drawingData), fps = 25)
# pendulumAnimation

animate(pendulumPlots, nframes = nrow(drawingData), fps = 25)
anim_save("doublePendulum.mp4")

system("ffmpeg -i doublePendulum.mp4 -vf 'setpts=1*PTS' doublePendulum_f.mp4")
