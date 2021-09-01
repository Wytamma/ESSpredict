# ESS predict 


```R
source('ESSpredict.R')
log <- read.logFile('logs/4.2.1-het-strict.aln.log')
target <- 200
burnin <- 0.1
ESSdf <- cumulativeESS(log, cores = 6, burnin = burnin)
plotESS(ESSdf)
```
![](figures/1.png)

```R
fit <- stepsUntil(ESSdf, target)
parameter_below_ESS <- setdiff(unique(ESSdf$parameter), unique(ESSdf[ESSdf$ESS > target,]$parameter))
for (stat in parameter_below_ESS) {
  p <- plotESS(
      ESSdf[ESSdf$parameter == stat,][ESSdf[ESSdf$parameter == stat,]$ESS < target, ], 
      fit=fit[fit$parameter == stat,])
  print(p)
}
```
![](figures/2.png)
