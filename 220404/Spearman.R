#!/usr/bin/env Rscript
library(ggplot2)
library(ggsci)
library(gridExtra)
library(foreach)
library(doParallel)
library(readr)

setwd("~/wyf/cancer_DNAm_distance")
Args <- commandArgs(T)
data_path <- Args[1]
output_path <- Args[2]

do_spearman <- function(Site, Exp) {
  test_out <- cor.test(
    rep((1:7), sum(Exp)),
    do.call(c, lapply(as.numeric(colnames(Exp[, Exp == 1])), function(x)
      as.numeric(Site[[x]][rownames(Exp), ]))),
    conf.level = 0.95,
    method = "spearman",
    exact = FALSE
  )
  TP <- test_out$p.value
  SIG = 0
  if (!is.na(TP)) {
    if (TP < 0.05) {
      SIG = 1
    }
  }
  TR <- as.numeric(test_out$estimate)
  SP <- array(dim = 10)
  SR <- array(dim = 10)
  SSIG <- NULL
  SUM_SIG = 0
  for (i in as.numeric(colnames(Exp[, Exp == 1]))) {
    test_out <- cor.test(
      c(1:7),
      as.numeric(Site[[i]][rownames(Exp), ]),
      conf.level = 0.95,
      method = "spearman",
      exact = FALSE
    )
    SP[i] <- test_out$p.value
    SR[i] <- as.numeric(test_out$estimate)
    if (!is.na(test_out$p.value)) {
      if (test_out$p.value < 0.1) {
        SUM_SIG = SUM_SIG + 1
        SSIG <- c(SSIG, 0.88)
      } else{
        SSIG <- c(SSIG, 0.28)
      }
    } else {
      SSIG <- c(SSIG, 0.28)
    }
  }
  if (SUM_SIG > 2 | SIG > 0) {
    draw_plot(Site, Exp, SSIG, TP, TR)
  }
  return(as.data.frame(cbind(
    rownames(Exp), sum(Exp), SUM_SIG, t(SP), t(SR), TP, TR
  )))
}

draw_plot <- function(Site, Exp, Sig, TP, TR) {
  figfile <-
    paste("png/Site", rownames(Exp), "_spearman.png", sep = "")
  p <- ggplot() +
    geom_point(aes(x = rep((1:7), sum(Exp)),
                   y = do.call(
                     c, lapply(as.numeric(colnames(Exp[, Exp == 1])), function(x)
                       as.numeric(Site[[x]][rownames(Exp),]))
                   ))) +
    geom_line(aes(
      x = rep((1:7), sum(Exp)),
      y = do.call(c, lapply(as.numeric(colnames(Exp[, Exp == 1])), function(x)
        as.numeric(Site[[x]][rownames(Exp),]))),
      size = as.factor(rep(sprintf(
        "%02d", as.numeric(colnames(Exp[, Exp == 1]))
      ), each = 7)),
      color = as.factor(rep(sprintf(
        "%02d", as.numeric(colnames(Exp[, Exp == 1]))
      ), each = 7))
    )) +
    scale_color_jco(name = "Patient") +
    scale_x_continuous(
      breaks = c(1:7),
      limits = c(0.8, 7.2),
      expand = c(0, 0),
      labels = c("T", "TE", "P5", "P10", "P15", "P20", "PN")
    ) +
    scale_size_manual(values = Sig,
                      guide = 'none') +
    xlab(label = "Location") +
    ylab(label = "Beta value") +
    ggtitle(label = paste("P-value=", TP, "\nrho=", TR, sep = "")) +
    theme(
      legend.background = element_blank(),
      legend.spacing = unit(0.1, units = "mm"),
      legend.key.size = unit(3.2, 0.2, units = "mm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "#000000"),
      axis.ticks = element_line(colour = "#000000"),
      axis.ticks.length.y = unit(2, units = "mm"),
      axis.ticks.length.x = unit(1, units = "mm"),
      axis.title = element_text(
        colour = "#000000",
        face = "bold",
        size = 14
      ),
      axis.text = element_text(
        colour = "#000000",
        face = "bold",
        size = 12
      )
    )
  png(figfile,
      width = 1200,
      height = 800,
      res = 200)
  print(p)
  dev.off()
}

write("==> Start read file!", stderr())

point_list <-
  read.csv(data_path,
           sep = "\t",
           header = F,
           row.names = 1)
colnames(point_list) <- (1:10)
point_list <-
  point_list[order(as.numeric(row.names(point_list))),]

P <- list()
for (i in (1:10)) {
  input_file <-
    paste("sitefilter15_P", sprintf("%02d", i), ".tsv", sep = "")
  dft <- read.csv(input_file,
                  sep = "\t",
                  header = F,
                  row.names = 1)
  colnames(dft) <- paste("P", (1:7), sep = "")
  P[[i]] <- data.frame(dft)
}

write("==> Start parallel!", stderr())
cl <- makeCluster(20)
registerDoParallel(cl)
temp <- foreach(
  ID = rownames(point_list),
  .packages = c("dplyr", "ggplot2", "readr", "gridExtra", "ggsci"),
  .inorder = F
) %dopar% {
  Site <- list()
  for (samp in 1:10) {
    Site[[samp]] <- round(P[[samp]][ID,], 0)
  }
  output <-
    do_spearman(Site, point_list[ID,])
  if (!is.null(output)) {
    write_delim(
      output,
      output_path,
      append = T,
      col_names = F,
      delim = "\t"
    )
  }
}
quit()