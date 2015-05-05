plot.fit.normalmixEM <- function(result, ...) {
	dmnorm = function(x, lambda, mu, sigma) {
		a = 0
		for (i in 1:length(lambda))
			a = a + lambda[i]*dnorm(x, mu[i], sigma[i])
		a
	}
	curve(dmnorm(x, result$lambda, result$mu, result$sigma), n=200, ...)
}