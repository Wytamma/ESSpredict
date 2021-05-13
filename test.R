source('ESSpredict.R')

log <- read.logFile('4.2.2-het-UCLN.aln.log')
essdf <- cumESS(log[1:16], cores = 4)
cutoff = 1e7
target <- 200
fit <- stepsUntil(essdf[essdf$state < cutoff,], target)
Statistic_below_ESS <- setdiff(unique(essdf$Statistic), unique(essdf[essdf$ESS > target,]$Statistic))
for (stat in unique(essdf$Statistic)) {
  p <- plotESS(essdf[essdf$Statistic == stat,][essdf[essdf$Statistic == stat,]$ESS < target, ], cutoff=cutoff, fit=fit[fit$Statistic == stat,])
  print(p)
}

l <- list()

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
