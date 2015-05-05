plot.boot.comp = function(result, B, sig, arbmean, arbvar) {
	if (length(result$p.values) != 2) {
		g = ceiling(sqrt(length(result$p.values)))
		par(mfrow = c(g, g))
	}
	else {
		par(mfrow = c(1, 2))
	}

	for (i in 1:length(result$p.values)) {
		chisq.df = ifelse(i == 1, 4, 3)
		if (!arbmean)  chisq.df = chisq.df - (i+1)
		if (!arbvar)   chisq.df = chisq.df - (i+1)
		xlim = c(0, max(result$log.lik[[i]], result$obs.log.lik[i]))
		ylim = c(0, max(max(hist(result$log.lik[[i]], plot=FALSE)$density),
					optimize(dchisq, xlim, df=chisq.df, maximum=TRUE)$objective))

		hist(result$log.lik[[i]], freq=FALSE, xlim=xlim, ylim=ylim, xlab='Bootstrap LRTS', main=paste(i, 'versus', i + 1, 'Components'))
		points(result$obs.log.lik[i], 0, col=4, pch=16)
		curve(dchisq(x, df=chisq.df), add=TRUE, n=1000, from=0, col=2)
		segments(qchisq(1-sig, df=chisq.df), 0, qchisq(1-sig, df=chisq.df), B, lty=2, col=2)
	}
}

boot.comp = function(y, x = NULL, N = NULL, max.comp = 2, B = 100, sig = 0.05, arbmean = TRUE, arbvar = TRUE, mix.type = c('logisregmix', 'multmix', 'mvnormalmix', 'normalmix', 'poisregmix', 'regmix', 'regmix.mixed', 'repnormmix'), hist = TRUE, ...) {

	(result = mixtools::boot.comp(y, x, N, max.comp, B, sig, arbmean, arbvar, mix.type, hist=FALSE, ...))

	# Let's create a more informative histogram
	if (mix.type == 'normalmix' && hist)
		plot.boot.comp(result, B, sig, arbmean, arbvar)

	result
}