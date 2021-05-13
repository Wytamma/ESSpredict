library(coda)
library(runner)
library(tidyverse)
library(ggplot2)
library(parallel)
library(lme4)

read.logFile <- function(path_to_log_file, index_col=NULL) {
  logFile <- read.csv(path_to_log_file, sep="\t", comment.char = '#')
  if (is.null(index_col)) {
    if ('state' %in% colnames(logFile)) {
      index_col = 'state'
    } else if ('Sample' %in% colnames(logFile)) {
      index_col = 'Sample'
    } else {
      error('index column not found in logFile, please set index_col')
    }
  }
  # rename
  names(logFile)[names(logFile) == index_col] <- "state"
  return(logFile)
}

cumESS <- function(logFile, burnin=0.1, sample=1, tidy=TRUE, cores=detectCores()) {
  # burn in
  if (burnin > 1 ) {
    logFile <- logFile[logFile[,"state"] >= burnin, ]
  } else {
  logFile <- logFile[logFile[,"state"] >= (tail(logFile[,"state"], n=1) * burnin), ]
  }
  # sample
  if (sample < 1) {
    logFile <- logFile[sample(nrow(logFile), sample*nrow(logFile), replace=FALSE), ]
    logFile <- logFile[order(as.numeric(rownames(logFile))), ]
  }

  cumESS_list <- mclapply(
    logFile[ , -which(names(logFile) %in% c("state")), drop=FALSE], # don't cal ESS on index
    runner::runner,
    f=function(x) if (length(x) == 1) 0 else coda::effectiveSize(x),
    mc.cores = cores) # doesn't work on windows unless cores = 1

  cumESSdf <- data.frame(cumESS_list)
  cumESSdf[,"state"] <- logFile[,"state"] #

  if (tidy) {
    cumESSdf <- cumESSdf %>%
      pivot_longer(!all_of("state"), names_to = "Statistic", values_to = "ESS")
  }
  return(cumESSdf)
}

plotESS <- function(tidyCumESSdf, fit=NULL, cutoff=NULL) {
  if (is.null(cutoff)) {
    p <- tidyCumESSdf %>%
      ggplot() + geom_point(aes(x = state, y = ESS, color = Statistic))
  } else {
    p <- ggplot() + geom_point(
      data=tidyCumESSdf %>%
        filter(state < cutoff),
      aes(x = state, y = ESS, color = Statistic))
    p <- p + geom_point(
      data=tidyCumESSdf %>%
        filter(state >= cutoff),
      aes(x = state, y = ESS), color = 'grey70')
  }
  if (!is.null(fit)) {
    p <- p +
      #geom_ribbon(data=fit, mapping = aes(x = fit, y=ESS, xmin = fit - se.fit, xmax = fit + se.fit, color = Statistic), fill = "grey70") +
      geom_line(data=fit, mapping = aes(x = fit, y=ESS, color = Statistic))
  }
  return(p)
}

stepsUntil <- function(tidyCumESSdf, ESStarget) {
  Statistics <- unique(tidyCumESSdf$Statistic)
  Statistic_above_ESS <- unique(tidyCumESSdf[tidyCumESSdf$ESS > ESStarget,]$Statistic)
  Statistic_below_ESS <- setdiff(Statistics, Statistic_above_ESS)

  above_ESS_df <- data.frame(list(statistic=Statistic_above_ESS, steps=rep.int(0, length(Statistic_above_ESS))))
  if (length(Statistic_below_ESS) == 0) {
    message('ESS targets reached!')
  }
  models <- lmList(state ~ ESS | Statistic, data=as.data.frame(tidyCumESSdf))
  fit <- predict(models, newdata = data.frame(ESS=0:ESStarget), se.fit = TRUE)
  fit$current.state <- tail(tidyCumESSdf$state, n=1)
  fit$ESS <- 0:ESStarget
  fit$estimated_steps <- fit$fit - fit$current.state

  return(fit)
}

