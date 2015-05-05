plot.mvnormalmix = function(x, y, lambda, mu, sigma, interactive=TRUE, lwd=0.1, ...) {
	if (interactive & !require(rgl))
		stop('You need to install the rgl package.')
	
	len = max(100*(max(x)-min(x)), 100*(max(y)-min(y)))

	x = seq(min(x), max(x), length.out=len)
	y = seq(min(y), max(y), length.out=len)
	z = outer(x, y, function(x, y) {
		m = cbind(x, y)
		n = rep(0, dim(m)[1])
		for (i in seq_along(lambda)) {
			n = n + dmvnorm(m, mu[[i]], sigma[[i]])*lambda[i]
		}
		return(n)
	})

	if (interactive) {
		height = (z - range(z)[1]) / diff(range(z))
		r.prop = height
		g.prop = 0
		b.prop = 1-height
		color  = rgb(r.prop, g.prop, b.prop, maxColorValue=1)

		persp3d(x, y, z, col=color, ...)
	} else {
		persp(x, y, z, lwd=lwd, ...)
	}
}

plot.mvnormalmixEM = function(result, interactive=TRUE, ...) {
	if (class(result) != 'mixEM' || result$ft != 'mvnormalmixEM')
		stop('Result invalid. Please use an output of mvnormalmixEM.')
	if (ncol(result$x) != 2)
		stop('Result invalid. Please use this function for **bivariate** Normals only.')

	plot.mvnormalmix(x=result$x[,1], y=result$x[,2], lambda=result$lambda, mu=result$mu, sigma=result$sigma, interactive=interactive, ...)
}