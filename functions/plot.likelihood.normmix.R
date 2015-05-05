plot.likelihood.normmix = function(data, muOne, muTwo, sigma, lambdas, n=100, columns=2, theta=30, phi=30, lwd=0.1, expand=0.5, ...) {
  par(mfrow=c(ceiling(length(lambdas)/columns), columns))

  x = seq(muOne[1], muOne[2], length.out=n)
  y = seq(muTwo[1], muTwo[2], length.out=n)

  for (lambda in lambdas) {
    par(mar=c(1,1,2,1))
    z = outer(x, y, function(x, y) {
      apply(cbind(x, y), 1, function(mu) {
        prod(dnorm(data, mu[1], sigma)*lambda + dnorm(data, mu[2], sigma)*(1-lambda))
      })
    })

    persp(x, y, z, xlab='mu_1', ylab='mu_2', zlab='L', main=paste('Pi_1 = ', lambda), theta=theta, phi=phi, lwd=lwd, expand=expand, ...)
  }
}

data  = as.numeric(faithful$eruptions)
muOne = c(min(data), max(data))
muTwo = muOne
sigma = sd(data)
lambdas = c(0.55,.0.75)
plot.likelihood.normmix(data, muOne, muTwo, sigma, lambdas, theta=90)