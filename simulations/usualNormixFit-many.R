source('../functions/usualNormixFit.R')
set.seed(2015)
i = 1
for (lambda in seq(0.1, 0.5, by=0.1)) {
	for (muDiff in seq(0.5, 3, by=0.33)) {
		for (sigmaDiff in seq(0, 2, by=0.33) {
			for (fitFn in c('nor1mix', 'mixtools')) {
				name = paste('mfit', fitFn, 'l', lambda, 'md', muDiff, 'sd', sigmaDiff, sep='_')
				if(fitFn == 'nor1mix') {
					fitFnArg = 'nor1mix'
				} else {
					fitFnArg = NULL
				}
				png(paste('images/', name, '.png', sep=''), width=1200*7/7, height=1200*4/7, pointsize=26)

				fit = usualNormixFit(lambda=c(lambda, 1-lambda), mu=c(0, muDiff), sigma=c(1, 1+sigmaDiff), fitFn = fitFnArg)

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