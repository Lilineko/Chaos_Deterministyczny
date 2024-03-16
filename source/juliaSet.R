library(viridis)
library(caTools)

# resolution
xPixels <- 240
yPixels <- 160

# range in a complex plane
xRange <- 2.60;
yRange <- 2.15;

xOffset = 0.75

# axes
x <- xRange * ((1 - xPixels / 2) : (xPixels / 2)) / xPixels - xOffset
y <- yRange * ((1 - yPixels / 2) : (yPixels / 2)) / yPixels

# color scaling function
scalingMandelbrot <- function(x) 1.0 - atan(5 * x**10)

#definition of the function calculating the Mandelbrot set
mandelbrot <- function(x, y, maxIteration = 100) {
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

divergence <- mandelbrot(x, y)
mandelbrotColorMap <- matrix(scalingMandelbrot(divergence), nrow = xPixels, byrow = FALSE)
image(mandelbrotColorMap, col = magma(100), axes = FALSE)

# now one can track a trajectory in the complex plane with respect to the Mandelbrot set
# this trajectory will produce a family of Julia sets that we would like to animate

# resolution
xJuliaPixels <- 800
yJuliaPixels <- 600

# range in a complex plane
xJuliaRange <- 4.0;
yJuliaRange <- 3.0;

# axes
x <- xJuliaRange * ((1 - xJuliaPixels / 2) : (xJuliaPixels / 2)) / xJuliaPixels
y <- yJuliaRange * ((1 - yJuliaPixels / 2) : (yJuliaPixels / 2)) / yJuliaPixels

# definition of the function calculating single Julia set
julia <- function(x, y, c, maxIteration = 100) {
  z0 <- outer(x, 1i * y, '+') %>% c()
  z0 <- array(z0, c(length(x), length(y)))
  result = array(0.0, c(length(x), length(y)))
  for(ix in 1:length(x)) {
    for(iy in 1:length(y)) {
      z = z0[ix, iy]
      iteration = 0
      while(iteration < maxIteration) {
        if(abs(z) > 2.0) {
          result[ix, iy] = 1.0 - (iteration / maxIteration)
          break
        }
        iteration = iteration  + 1
        z = z * z + c
      }
    }
  }
  return(result)
}

getCircle <- function(x, r, x0 = 0, y0 = 0) y0 - sqrt(r**2 - (x-x0)**2)

# setting a path in the complex plane
c <- seq(-1.25, -1, length.out = 30)
c = c(c, seq(-1.0, -0.5, length.out = 60) + 1i * seq(0.0, 0.666, length.out = 60))
c = c(c, seq(-0.5, 0.1, length.out = 60) + 1i * seq(0.666, 0.666, length.out = 60))
c = c(c, seq(0.1, 0.4, length.out = 40) + 1i * seq(0.666, 0.366, length.out = 40))
c = c(c, seq(0.4, 0.4, length.out = 90) + 1i * seq(0.366, -0.45, length.out = 90))
c = c(c, seq(0.4, 0.0, length.out = 70) + 1i * seq(-0.45, -0.8, length.out = 70))
tempx = seq(0.0, -0.75, length.out = 140)
tempy = getCircle(tempx, seq(0.75, 0.75, length.out = length(tempx)), x0 = 0, y0 = -0.05)
c = c(c, tempx + 1i * tempy)
tempx = seq(-0.75, -1.25, length.out = 80)
tempy = getCircle(tempx, seq(0.25, 0.25, length.out = length(tempx)), x0 = -1.0, y0 = -0.05)
c = c(c, tempx + 1i * tempy)
c = c(c, seq(-1.25, -1.25, length.out = 10) + 1i * seq(-0.05, 0.0, length.out = 10))

scalingJulia <- function(x) (1.0 - atan(10 * x**4) + 0.4) / 1.4

curveColor = 0.5
insertion = mandelbrotColorMap
xShift = xJuliaPixels - xPixels
yShift = yJuliaPixels - yPixels
juliaSets <- array(0, c(length(x), length(y), length(c)))
for(it in 1:length(c)) {
  set = abs(julia(x, y, c[it]))
  juliaSets[ , , it] = scalingJulia((set - min(set)) / (max(set) - min(set)))

  # inserting the Mandelbrot set with the path
  im = Re(c[it])
  re = Im(c[it])
  xPos = round(xPixels * (1 / 2 + (im + xOffset) / xRange))
  yPos = round(yPixels * (1 / 2 + re / yRange))
  insertion[(xPos-1):(xPos+1), (yPos-1):(yPos+1)] = curveColor
  for(ix in 1:xPixels) {
    for(iy in 1:yPixels) {
      juliaSets[ix + xShift, iy + yShift, it] = max(c(juliaSets[ix + xShift, iy + yShift, it], insertion[ix, iy]))
    }
  }
}
write.gif(juliaSets, "julia.gif", col = magma(256), flip = TRUE, scale = "smart")

