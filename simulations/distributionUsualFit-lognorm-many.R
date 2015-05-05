source('../functions/distributionUsualFit.R')
i = 1
for (meanlog in c(100, 1000, 5000, 10000, 15000, 20000, 30000, 50000)) {
	for (sdlog in c(0.25, 0.4, 1, 2, 10)) {
		for (k in 2:5) {
			name = paste('lognorm_fit', 'k', k, 'meanlog', meanlog, 'sdlog', sdlog, sep='_')
			png(paste('images/', name, '.png', sep=''), units='px', width=1200*7/7, height=1200*4/7, pointsize=26)

			set.seed(2015)
			fit = distributionUsualFit(s=200, k=k, xlim=c(0, 2*meanlog), meanlog=log(meanlog), sdlog=sdlog)

			dev.off()
			assign(name, fit)
			save(list=name, name, file=paste('fits/', name, '.RData', sep=''))
			print(paste('Done number', i))
			print(paste('Done', name))
			i <<- i + 1
		}
	}
}