distributionUsualFit = function(n=1000, s=1000, k=NULL, theorDistr='lognorm', randomFn=NULL, densityFn=NULL, xlim=NULL, fit=NULL, fitFn=NULL, usualFn=NULL, ...) {
	# Require packages
	if (!require(nor1mix) || !require(parallel))
		stop('You need nor1mix and parallel packages installed.')
	# Number of components
	if (is.null(k) || is.null(xlim) || length(xlim) != 2)
		stop('Please pass valid k and xlim.')
	# Construct theoretical density
	if (!is.function(densityFn) || !is.function(randomFn)) {
		if (theorDistr == 'lognorm') {
			densityFn = dlnorm
			randomFn = rlnorm
		} else {
			stop('Please pass valid theorDistr or randomFn and densityFn.')
		}
	}
	if (is.null(fit)) {
		# Allow custom fit function - for package comparisons
		if (!is.function(fitFn)) {
			if (is.null(fitFn) || fitFn == 'nor1mix') {
				fitFn = function(x, k, ...) {
					a = norMixEM(x, m=k, maxiter=10000)
					return(list(mu = a[,'mu'], sigma = a[,'sigma'], lambda = a[,'w']))
				}
			} else if (fitFn == 'mixtools') {
				fitFn = function(x, k, ...) normalmixEM(x, k=k)[2:4]
			} else {
				stop('Please pass a valid fit function.')
			}
		}
		fit = mclapply(1:s, function(i, ...){
			d = randomFn(n, ...)
			a = fitFn(d, k)
			# Order the components by mu 
			lapply(a, function(x) x[order(a$mu)] )
		}, ...)
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
	if (any(is.na(uf$sigma)) || any(uf$sigma <= 0)) {
		print('Invalid sigma, skipping.')
		return(FALSE)
	}
	# Construct "usual fit" mixture
	fnm = norMix(mu=uf$mu, sigma=uf$sigma, w=uf$lambda)
	dnm = list(fit = dnorMixL(fnm))
	# Get X and Y coordinates
	xTheoretical = seq(xlim[1], xlim[2], length.out=1000)
	yTheoretical = densityFn(xTheoretical, ...)
	# Plot real density and "usual fit"
	par(mar=c(4,2,2,1))
	ylim = c(min(yTheoretical, dnm$fit$y), max(yTheoretical, dnm$fit$y))
	plot(yTheoretical ~ xTheoretical, xlim = xlim, ylim = ylim, main = 'Mixture density', xlab = NA, ylab = NA, type='l')
	lines(fnm, col=2)
	for(j in 1:k)
		curve(fnm[j,'w']*dnorm(x, fnm[j,'mu'], fnm[j,'sigma']), col='#ff8080', lty=2, add=TRUE)
	legend('topright', legend=c('real p.d.f.', 'usual fit', 'fit comp.'), col=c(1, 2, '#ff8080'), lwd=1, lty=c(1,1,2))
	
	# Return a list of everything
	return(list(n=n, s=s, k=k, theorDistr=theorDistr, randomFn=randomFn, densityFn=densityFn, fit=fit, fitFn=fitFn, usualFn=usualFn, uf=uf, est=est))
}