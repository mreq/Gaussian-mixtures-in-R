plot.likelihood.norm <- function(data, meanRange, sdRange, interactive=TRUE, n=100, ...) {
  if (interactive & !require(rgl))
    return(print('You need to install the rgl package.'))

  n <- length(data)

  x <- seq(meanRange[1], meanRange[2], length.out=n)
  y <- seq(sdRange[1],   sdRange[2],   length.out=n)
  z <- outer(x, y, function(x, y) {
    apply(cbind(x, y), 1, function(theta) {
      prod((dnorm(data, theta[1], theta[2])))
    })
  })
  z <- z/max(z)

  r.prop <- z
  g.prop <- 0
  b.prop <- 1-z
  color  <- rgb(r.prop, g.prop, b.prop, maxColorValue=1)

  if (interactive) {
    persp3d(x, y, z, col=color, xlab='mu', ylab='sigma', zlab='Relative L', ...)
  } else {
    persp(x, y, z, xlab='mu', ylab='sigma', border=NA, shade=0.6, col=2, ...)
  }
}

data      <- as.numeric(ldeaths)
meanRange <- c(min(data), max(data))
sdRange   <- c(0.5*sd(data), 2*sd(data))
plot.likelihood.norm(data, meanRange, sdRange)