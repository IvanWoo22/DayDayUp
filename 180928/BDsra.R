#!/usr/bin/r

args <- commandArgs(T)


downlist <- read.table(args[1],header = T)
setwd(args[2])


for (i in downlist$sra_ftp){
  url = paste("aria2c -x 4 -c ftp://",i,sep = "",collapse="")
  system(url)
}
