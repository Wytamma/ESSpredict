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
