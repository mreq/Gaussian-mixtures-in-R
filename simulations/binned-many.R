source('../functions/usualNormixFit.R')
library(bda)
fitFn = function(x, k, ...) {
	# Calculate binning breaks
	bn = binning(x, breaks=c(-Inf, hist(x, plot=FALSE)$breaks, Inf))
	a = fit.mixnorm(bn, k)
	return(list(lambda=a$p, mu=a$mu, sigma=a$s))
}
set.seed(2015)
i = 1
for (lambda in seq(0.1, 0.5, by=0.1)) {
	for (muDiff in seq(0.5, 3, by=0.33)) {
		for (sigmaDiff in seq(0, 2, by=0.33)) {
			name = paste('binfit', 'l', lambda, 'md', muDiff, 'sd', sigmaDiff, sep='_')
			png(paste('images/', name, '.png', sep=''), width=1200*7/7, height=1200*4/7, pointsize=26)

			fit = usualNormixFit(k=2, lambda=c(lambda, 1-lambda), mu=c(0, muDiff), sigma=c(1, 1+sigmaDiff), fitFn=fitFn)

			dev.off()
			assign(name, fit)
			save(list=name, name, file=paste('fits/', name, '.RData', sep=''))
			rm(list=paste(name))
			print(paste('Done number', i))
			i <<- i + 1
			print(paste('Done', name))
		}
	}
}