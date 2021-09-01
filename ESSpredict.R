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

cumulativeESS <- function(logFile, burnin=0.1, sample=1, tidy=TRUE, cores=detectCores() - 1, k = 1) {
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

  ESS <- function(x) {
    if (length(x) %% k == 0 && length(x) != 1) {
      return(coda::effectiveSize(x))
    } else {
      return(NA)
    }
  }

  cumulativeESS_list <- mclapply(
    logFile[ , -which(names(logFile) %in% c("state")), drop=FALSE], # don't cal ESS on index
    runner::runner,
    f=ESS,
    mc.cores = cores) # doesn't work on windows unless cores = 1

  cumulativeESSdf <- data.frame(cumulativeESS_list)
  cumulativeESSdf[,"state"] <- logFile[,"state"] #

  if (tidy) {
    cumulativeESSdf <- cumulativeESSdf %>%
      pivot_longer(!all_of("state"), names_to = "parameter", values_to = "ESS")
  }
  return(cumulativeESSdf)
}

plotESS <- function(tidycumulativeESSdf, fit=NULL, cutoff=NULL) {
  if (is.null(cutoff)) {
    p <- tidycumulativeESSdf %>%
      drop_na() %>%
      ggplot() + geom_point(aes(x = state, y = ESS, color = parameter))
  } else {
    p <- ggplot() + geom_point(
      data=tidycumulativeESSdf %>%
        filter(state < cutoff),
      aes(x = state, y = ESS, color = parameter))
    p <- p + geom_point(
      data=tidycumulativeESSdf %>%
        filter(state >= cutoff),
      aes(x = state, y = ESS), color = 'grey70')
  }
  if (!is.null(fit)) {
    p <- p +
      geom_line(data=fit, mapping = aes(x = fit, y=ESS, color = parameter))
  }
  return(p)
}

stepsUntil <- function(tidycumulativeESSdf, ESStarget) {
  parameters <- unique(tidycumulativeESSdf$parameter)
  parameter_above_ESS <- unique(tidycumulativeESSdf[tidycumulativeESSdf$ESS > ESStarget,]$parameter)
  parameter_below_ESS <- setdiff(parameters, parameter_above_ESS)

  above_ESS_df <- data.frame(list(parameter=parameter_above_ESS, steps=rep.int(0, length(parameter_above_ESS))))
  if (length(parameter_below_ESS) == 0) {
    message('ESS targets reached!')
  }
  models <- lmList(state ~ ESS | parameter, data=as.data.frame(tidycumulativeESSdf))
  fit <- predict(models, newdata = data.frame(ESS=0:ESStarget), se.fit = TRUE)
  fit$current.state <- tail(tidycumulativeESSdf$state, n=1)
  fit$ESS <- 0:ESStarget
  fit$estimated_steps <- fit$fit - fit$current.state

  return(fit)
}

parametersBelow <- function(ESSdf, target) {
  return(setdiff(unique(ESSdf$parameter), unique(ESSdf[ESSdf$ESS > target,]$parameter)))
}

