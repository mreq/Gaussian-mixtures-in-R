source('../functions/usualNormixFit.R')
i = 1
for (lambda in seq(0.1, 0.5, by=0.1)) {
	for (muDiff in seq(0.5, 3, by=0.5)) {
		for (sigmaDiff in seq(0, 2, by=0.5)) {
			for (fitFn in c('nor1mix', 'mixtools')) {
				name = paste('overfit', fitFn, 'l', lambda, 'md', muDiff, 'sd', sigmaDiff, sep='_')
				if(fitFn == 'nor1mix') {
					fitFnArg = 'nor1mix'
				} else {
					fitFnArg = NULL
				}
				png(paste('images/', name, '.png', sep=''), units='px', width=1200*7/7, height=1200*4/7, pointsize=26)
				# Cairo(paste('images/', name, '.png', sep=''), type='png', units='px', width=1200*7/7, height=1200*4/7, pointsize=26, antialias=TRUE)
        
				set.seed(2015)
				fit = usualNormixFit(k=3, lambda=c(lambda, 1-lambda), mu=c(0, muDiff), sigma=c(1, 1+sigmaDiff), fitFn = fitFnArg)

				dev.off()
				assign(name, fit)
				save(list=name, name, file=paste('fits/', name, '.RData', sep=''))
				print(paste('Done number', i))
				i <<- i + 1
				print(paste('Done', name))
			}
		}
	}
}