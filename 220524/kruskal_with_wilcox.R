#!/usr/bin/env Rscript
library(tidyverse)
library(rstatix)
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
    tibble(
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
    res.kru1 <- kruskal.test(beta ~ loc, data = selfesteem)$p.value
    tmp <-
      cbind(rawdat[i, 1], rawdat[i, 2], res.kru1)
    
  }, error = function(e) {
    tmp[1, 1] <<- rawdat[i, 1]
    tmp[1, 2] <<- rawdat[i, 2]
    tmp[1, 3] <<- NA
  })
  tryCatch({
    res.kru2 <- kruskal.test(beta ~ patient, data = selfesteem)$p.value
    tmp <-
      cbind(tmp, res.kru2)
  }, error = function(e) {
    tmp[1, 4] <<- NA
  })
  tmp <- as.data.frame(tmp)
  pwc1 <-
    selfesteem %>% wilcox_test(
      beta ~ loc,
      p.adjust.method = "none",
      paired = T,
      alternative = "greater"
    )
  tmp1 <- as.data.frame(t(pwc1$p.adj))
  colnames(tmp1) <- paste(pwc1$group1,pwc1$group2,"greater",sep = "_")
  pwc2 <-
    selfesteem %>% wilcox_test(
      beta ~ loc,
      p.adjust.method = "none",
      paired = T,
      alternative = "less"
    )
  tmp2 <- as.data.frame(t(pwc2$p.adj))
  colnames(tmp2) <- paste(pwc2$group1,pwc2$group2,"less",sep = "_")
  tmp <- cbind(tmp,tmp1,tmp2)
  write_delim(
    tmp,
    output_path,
    append = T,
    col_names = F,
    delim = "\t"
  )
}
quit()