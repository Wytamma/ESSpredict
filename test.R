source('ESSpredict.R')

log <- read.logFile('/Users/wytamma/programming/typhi-in-fiji/skygrid/4.2.2-skygrid.aln.log')
cutoff <- 1e8
target <- 200
burnin <- 0.1
essdf <- cumESS(log, cores = 6, burnin = burnin)
plotESS(essdf)
fit <- stepsUntil(essdf, target)
Statistic_below_ESS <- setdiff(unique(essdf$Statistic), unique(essdf[essdf$ESS > target,]$Statistic))
for (stat in Statistic_below_ESS) {
  p <- plotESS(essdf[essdf$Statistic == stat,][essdf[essdf$Statistic == stat,]$ESS < target, ], cutoff=cutoff, fit=fit[fit$Statistic == stat,])
  print(p)
}

l <- list()

# do a first passs ESS of the full chain to see which statistics have EEE < target

for (r in seq(0, 1, 0.1)){
  xt <- 0
  for(i in 2:1000)  xt[i] <- r * xt[i-1] + rnorm(1, 0, 10)
  l[[paste0('r-', r)]] <- xt
}

l[['state']] <- 1:1000
autodf <- cumESS(as.data.frame(l), cores = 4)


for (r in unique(autodf$Statistic)) {
  x = autodf[autodf$Statistic == r,]$state
  y = autodf[autodf$Statistic == r,]$ESS
  print(paste0(r,' = ', (cor(x, y)[1] ^2 )))
}
