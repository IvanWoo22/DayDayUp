#!/usr/bin/env Rscript
library(tidyverse)
library(rstatix)
library(dplyr)
library(ggpubr)

Args <- commandArgs(T)
data_path <- Args[1]
output_path <- Args[2]

rawdat <-
  read.csv(data_path,
           header = F,
           sep = "\t",
           check.names = F)

for (i in 1:nrow(rawdat)) {
  selfesteem <-
    data.frame(
      beta = as.numeric(rawdat[i, 3:ncol(rawdat)]),
      loc = factor(
        rep(c(
          "TC", "TE", "P5", "P10", "P15", "P20", "PN"
        ), 12),
        levels = c("TC", "TE", "P5", "P10", "P15", "P20", "PN")
      ),
      patient = rep(paste("PA", sprintf("%02d", c(
        1, 5:15
      )), sep = ""), each = 7)
    )
  tryCatch({
    res.aov1 <-
      summary(aov(beta ~ loc, data = selfesteem))[[1]]$`Pr(>F)`[1]
    tmp <-
      cbind(rawdat[i, 1], rawdat[i, 2], res.aov1[1])
  }, error = function(e) {
    tmp[1, 1] <<- rawdat[i, 1]
    tmp[1, 2] <<- rawdat[i, 2]
    tmp[1, 3] <<- NA
  })
  tryCatch({
    res.aov2 <-
      summary(aov(beta ~ patient, data = selfesteem))[[1]]$`Pr(>F)`[1]
    tmp <-
      cbind(tmp, res.aov2[1])
  }, error = function(e) {
    tmp[1, 4] <<- NA
  })
  tmp <- as.data.frame(tmp)
  write_delim(
    tmp,
    output_path,
    append = T,
    col_names = F,
    delim = "\t"
  )
}
quit()