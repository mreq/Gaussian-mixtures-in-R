usualNormixFit = function(n=1000, s=1000, lambda=c(0.5, 0.5), mu=c(0, 1.75), sigma=c(1, 1), k=NULL, fit=NULL, fitFn=NULL, usualFn=NULL, smallPlots=TRUE) {
	# Validate arguments
	if (length(lambda) < 2 || length(lambda) != length(mu) || length(lambda) != length(sigma))
		stop('Please pass valid arguments lambda, mu and sigma.')
	# Require packages
	if (!require(mixtools) || !require(nor1mix) || !require(parallel))
		stop('You need mixtools, nor1mix, and parallel packages installed.')
	# Make lambda sum to one
	if (sum(lambda) != 1)
		lambda = lambda/sum(lambda)
	# Number of components
	if (is.null(k))
		k = length(lambda)
	# Construct theoretical mixture
	nm = norMix(mu=mu, sigma=sigma, w=lambda)
	# Allow saved fit - skip fitting and produce a plot only if fit is present
	if (is.null(fit)) {
		# Allow custom fit function - for package comparisons
		# defaults to mixtools::normalmixEM
		if (!is.function(fitFn)) {
			if (is.null(fitFn) || fitFn == 'mixtools') {
				fitFn = function(x, k, ...) normalmixEM(x, k=k)[2:4]
			} else if (fitFn == 'nor1mix') {
				fitFn = function(x, k, ...) {
					a = norMixEM(x, m=k, maxiter=10000)
					if (attr(a, "converged")) {
						return(list(mu = a[,'mu'], sigma = a[,'sigma'], lambda = a[,'w']))
					} else {
						return(list(mu = rep(NA, k), sigma = rep(NA, k), lambda = rep(NA, k)))
					}
				}
			} else {
				stop('Please pass a valid fit function.')
			}
		}
		fit = mclapply(1:s, function(x){
			rn = rnorMix(n, nm)
			a = fitFn(rn, k, mu, sigma, lambda)
			# Order the components by mu 
			lapply(a, function(x) x[order(a$mu)] )
		}, mc.silent=TRUE)
	}
	# Get estimates
	est = list(
		lambda = t(sapply(fit, function(a) a$lambda)),
		mu =     t(sapply(fit, function(a) a$mu)),
		sigma =  t(sapply(fit, function(a) a$sigma))
	)
	# Allow custom "usual fit" function
	# defaults to trimmed mean
	if (!is.function(usualFn)) {
		if (is.null(usualFn) || usualFn == 'mean') {
			usualFn = function(x)  mean(x, trim=0.2)
		} else if (usualFn == 'mode') {
			usualFn = function(x){
				y = round(x, 2)
				as.numeric(names(table(y))[which.max(table(y))])
			}
		}
	}
	# Usual fit
	uf = list(
		lambda = apply(est$lambda, 2, usualFn),
		mu =     apply(est$mu, 2, usualFn),
		sigma =  apply(est$sigma, 2, usualFn)
	)
	# Construct "usual fit" mixture
	fnm = norMix(mu=uf$mu, sigma=uf$sigma, w=uf$lambda)
	if (smallPlots) {
		# How many small graphs will be produced?
		#  = number of parameters
		smallPlots = sum(sapply(uf, length))
		# Big plot spans across 4 rows
		layoutMatrix = matrix(c(1,1,1,1,2,2,3,4,2,2,5,6,2,2,7,8), ncol=4, byrow=TRUE)
		while ((max(layoutMatrix)-2) < smallPlots) {
			i = max(layoutMatrix)
			layoutMatrix = rbind(layoutMatrix, c(i+1, i+2, i+3, i+4))
		}
		# Small height for heading
		layoutHeights = c(1, rep(3, nrow(layoutMatrix)-1))
	} else {
		layoutMatrix = matrix(c(1,2), ncol=1)
		# Small height for heading
		layoutHeights = c(1, 6)
	}
	layout(layoutMatrix, heights=layoutHeights)
	# Plot heading
	par(mar=c(0,0,0,0))
	# Create heading text
	headingText = c('"s = ', s, ' ; n = ', n, ' ; "')
	for (i in seq_along(lambda))
		headingText = c(headingText, paste('~lambda[', i, ']~"="~"', round(lambda[i], 2), '"~"; "', sep=''))
	for (i in seq_along(mu))
		headingText = c(headingText, paste('~mu[', i, ']~"="~"', round(mu[i], 2), '"~"; "', sep=''))
	for (i in seq_along(sigma))
		headingText = c(headingText, paste('~sigma[', i, ']~"="~"', round(sigma[i], 2), ifelse(i == length(sigma), '"', '"~"; "'), sep=''))
	headingTextArg = c(headingText, '')
	names(headingTextArg)[length(headingTextArg)] = 'sep'
	headingTextArg = as.list(headingTextArg)
	(headingText = do.call(paste, headingTextArg))
	plot.new()
	text(0.5, 0.5, parse(text=headingText), cex=1.5, font=2)
	# Get X and Y coordinates
	dnm = list(real = dnorMixL(nm), fit = dnorMixL(fnm))
	# Plot real density and "usual fit"
	par(mar=c(4,2,2,1))
	xlim = c(min(dnm$real$x, dnm$fit$x), max(dnm$real$x, dnm$fit$x))
	ylim = c(min(dnm$real$y, dnm$fit$y), max(dnm$real$y, dnm$fit$y))
	plot(dnm$real, xlim = xlim, ylim = ylim, main = 'Mixture density', xlab = NA, ylab = NA, type='l')
	lines(fnm, col=2)
	for(j in 1:length(lambda))
		curve(nm[j,'w']*dnorm(x, nm[j,'mu'], nm[j,'sigma']), col='#808080', lty=2, add=TRUE)
	for(j in 1:k)
		curve(fnm[j,'w']*dnorm(x, fnm[j,'mu'], fnm[j,'sigma']), col='#ff8080', lty=2, add=TRUE)
	legend('topright', legend=c('real p.d.f.', 'real comp.', 'usual fit', 'fit comp.'), col=c(1,'#808080', 2, '#ff8080'), lwd=1, lty=c(1,2,1,2))
	if (smallPlots) {
		# Create small plots
		for (param in names(est)) {
			d = est[[param]]
			if (param == 'lambda') {
				breaks = seq(0, 1, length.out=30)
			} else {
				breaks = seq(min(d)-0.01, max(d)+0.01, length.out=30)
			}
			for (i in 1:k) {
				par(mar=c(4,2,2,1))
				xlab = paste(param, '[', i, ']', sep='')
				hist(d[,i], breaks=breaks, xlab=NULL, main=parse(text=xlab))
				segments(x0=eval(parse(text=xlab)), y0=0, y1=n, lty=2, col=3)
				segments(x0=uf[[param]][i], y0=0, y1=n, lty=2, col=4)
			}
		}
	}
	# Return a list of everything
	return(list(n=n, s=s, lambda=lambda, mu=mu, sigma=sigma, fit=fit, fitFn=fitFn, usualFn=usualFn, nm=nm, fnm=fnm, est=est))
}