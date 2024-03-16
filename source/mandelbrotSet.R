library(tidyverse)
library(ggraph)
library(viridis)

# resolution
xPixels <- 1200
yPixels <- 800

# range in a complex plane
xRange <- 3.6;
yRange <- 2.4;

# axes
x <- xRange * ((1 - xPixels / 2) : (xPixels / 2)) / xPixels - 0.5
y <- yRange * ((1 - yPixels / 2) : (yPixels / 2)) / yPixels

#definition of the function calculating the Mandelbrot set
mandelbrot <- function(x, y, maxIteration) {
  c = outer(x, 1i * y, '+') %>% c()        # creating a complex grid c
  result = rep(0.0, length(c))             # container for resulting data
  
  # for each point in the complex plane we check if it belongs to the Mandelbrot set
  for(ic in 1:length(c)) {
    z = c[ic]
    iteration = 0
    while(iteration < maxIteration) {
      if(abs(z) > 2.0) {
        result[ic] = 1.0 - (iteration / maxIteration)
        break
      }
      iteration = iteration  + 1
      z = z * z + c[ic]
    }
  }
  # returned value is a divergance speed that one may use as a color function
  # value 0 means that the sequence converges, i.e. point belongs to the Mandelbrot set
  return(result)
}

divergence <- mandelbrot(x, y, 100)

colorMap <- matrix(atan(5 * divergence**10), nrow = xPixels, byrow = FALSE)
image(colorMap, col = magma(100), axes = FALSE)
