#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(getopt)
  library(tidyverse)
  library(survival)
  library(survminer)
  library(timeROC)
  library(ggplot2)
  library(scales)
  library(grid)
  library(gridExtra)
  library(extrafont)
})

spec = matrix(
  c(
    "help",
    "h",
    0,
    "logical",
    "brief help message",

    "simple",
    "s",
    0,
    "logical",
    "formula as the title",

    "infile",
    "i",
    1,
    "character",
    "input filename",

    "formula",
    "f",
    1,
    "character",
    "formula filename",

    "outdir",
    "o",
    1,
    "character",
    "output filename",

    "threshold",
    "t",
    1,
    "integer",
    "threshold of time for censoring"
  ),
  byrow = TRUE,
  ncol = 5
)
opt = getopt(spec)

if (!is.null(opt$help)) {
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}

if (is.null(opt$infile)) {
  cat("--infile is need\n")
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}

if (is.null(opt$formula)) {
  cat("--formula is need\n")
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}

if (is.null(opt$threshold)) {
  opt$threshold <- 1825
}

if (is.null(opt$outdir)) {
  opt$outdir <- "."
}

# col_type = cols() suppress the output
tbl <- read_tsv(opt$infile, col_type = cols())

# marker, coef, prev_p, median, hr
formula <- read_tsv(opt$formula, col_type = cols())
formula <- formula[, c(1, 2, 3, 4, 5)]

# status: 1 for event and 0 for censored
# right censored
tbl$status <- ifelse(tbl$time >= opt$threshold, 0, tbl$status)

# functions
plot_roc <- function(tbl, predictor) {
  tbl_roc <-
    data.frame(
      time = tbl$time,
      status = tbl$status,
      predictor,
      null = numeric(nrow(tbl))
    )
  tbl_roc <-
    tbl_roc[which(!is.na(tbl_roc$time) &
                    !is.na(tbl_roc$status) &
                    !is.na(tbl_roc$predictor)),]

  roc <- timeROC(
    T = tbl_roc$time,
    delta = tbl_roc$status,
    marker = tbl_roc$predictor,
    cause = 1,
    weighting = "marginal",
    times = opt$threshold,
    iid = TRUE
  )

  roc_null <- timeROC(
    T = tbl_roc$time,
    delta = tbl_roc$status,
    marker = tbl_roc$null,
    cause = 1,
    weighting = "marginal",
    times = opt$threshold,
    iid = TRUE
  )

  # timeROC add a time of 0 to the first
  # named num, use [[ ]] to extract value
  rocauc <- signif(roc$AUC[[2]], digits = 3)

  rocp <- compare(roc, roc_null)$p_values_AUC[[2]]
  p_value <- signif(rocp, digits = 2)

  aucci <- confint(roc, level = 0.95)$CI_AUC
  cilow <- signif(aucci[[1]] / 100.0, digits = 2)
  cihigh <- signif(aucci[[2]] / 100.0, digits = 2)

  # plot setup
  breaks = seq(0, 1, 0.2)
  plot <- ggplot() +
    geom_segment(aes(
      x = 0,
      y = 0,
      xend = 1,
      yend = 1
    )) +
    scale_x_continuous(
      name = "1 - Specificity",
      limits = c(0, 1),
      breaks = breaks,
      expand = c(0.02, 0.02)
    ) +
    scale_y_continuous(
      name = "Sensitivity",
      limits = c(0, 1),
      breaks = breaks,
      expand = c(0.02, 0.02)
    ) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())

  plot <- plot +
    geom_text(
      data = data.frame(),
      aes(
        label = paste(
          "P = ",
          p_value,
          "\nAUC = ",
          rocauc,
          "\n95% CI: ",
          cilow,
          "-",
          cihigh,
          sep = ""
        ),
        x = 1,
        y = 0
      ),
      hjust = 1,
      vjust = 0,
      size = 3
    )

  # steps
  sens_spec <- data.frame(spec = roc$FP[, 2], sens = roc$TP[, 2])
  plot <- plot +
    geom_step(aes(x = spec, y = sens), data = sens_spec, lwd = 0.8)

  return(plot)
}

plot_km <- function(tbl, fit, diff) {
  kmp <- pchisq(diff$chisq, length(diff$n) - 1, lower.tail = F)
  kmp <- signif(kmp, digits = 2)

  hr <- (diff$obs[1] / diff$exp[1]) / (diff$obs[2] / diff$exp[2])
  cilow = exp(log(hr) - qnorm(0.975) * sqrt(1 / diff$exp[1] + 1 / diff$exp[2]))
  cihigh = exp(log(hr) + qnorm(0.975) * sqrt(1 / diff$exp[1] + 1 / diff$exp[2]))
  hr <- signif(hr, digits = 2)
  cilow <- signif(cilow, digits = 2)
  cihigh <- signif(cihigh, digits = 2)

  ggsurv <- ggsurvplot(
    fit,
    data = tbl,
    legend.labs = c("High-risk", "Low-risk"),
    linetype = c("solid", "solid"),
    size = 0.8,
    # line width
    palette = c("black", "grey50"),
    censor.shape = "|",
    censor.size = 2,
    xlab = "Time (days)",
    ylab = "Proportion Surviving",
    pval = FALSE,
    ggtheme = theme_bw()
  )
  plot <- ggsurv$plot

  plot <- plot +
    geom_text(
      data = data.frame(),
      aes(
        label = paste("P = ", kmp, "\nHR = ", hr, "\n95% CI: ", cilow, "-", cihigh, sep =
          ""),
        x = 0,
        y = 0
      ),
      hjust = 0,
      vjust = 0,
      size = 3
    )

  breaks = seq(0, 1, 0.2)
  plot <- plot +
    theme(aspect.ratio = 1) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = breaks,
      expand = c(0.02, 0.02)
    ) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(legend.background = element_rect(colour = "transparent", fill = "transparent")) +
    theme(
      legend.justification = c(1, 0),
      legend.position = c(1, 0),
      legend.direction = "vertical"
    ) +
    theme(legend.title = element_blank())

  return(plot)
}

do_plot <- function(tbl, row) {
  #----------------------------#
  # get parameters
  #----------------------------#
  markers <- strsplit(as.character(row[1]), "+", fixed = T)[[1]]
  coefs <- strsplit(as.character(row[2]), ",", fixed = T)[[1]]
  coefs <- as.numeric(coefs)
  score_median <- as.numeric(row[4])

  figfile <- as.character(row[1])
  figfile <- str_replace_all(figfile, "\\+", "_")
  figfile <- str_interp("${opt$outdir}/${figfile}.pdf")

  write(figfile, stderr())

  #----------------------------#
  # K-M
  #----------------------------#
  predictor <- numeric(nrow(tbl))
  for (i in 1:length(markers)) {
    predictor <- predictor + coefs[i] * tbl[[markers[i]]]
  }

  tbl$group <- ifelse(predictor < score_median, 1, 0)
  if (all(table(tbl$group) > 5) &&
    length(table(tbl$group)) > 1) {
    diff <- tryCatch({
      survdiff(Surv(time, status) ~ group,
               data = tbl,
               rho = 0)
    },
      error = function(cond) {
        message(cond, "\n")
        return(NA)
      })

    fit <- tryCatch({
      survfit(Surv(time, status) ~ group, data = tbl)
    },
      error = function(cond) {
        message(cond, "\n")
        return(NA)
      })

    #----------------------------#
    # plots
    #----------------------------#
    pl_km <- plot_km(tbl, fit, diff)
    pl_roc <- plot_roc(tbl, predictor)
  } else
  {
    pl_km <- ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = "Insufficient group sizes",
        size = 5,
        hjust = 0.5
      ) +
      theme_void()

    pl_roc <- ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = "Insufficient group sizes",
        size = 5,
        hjust = 0.5
      ) +
      theme_void()
  }
  pdf(
    file = figfile,
    family = "Helvetica",
    width = 7,
    height = 4,
    useDingbats = FALSE
  )

  if (is.null(opt$threshold)) {
    opt$threshold <- 1825
  }
  title <- if (is.null(opt$simple)) {
    paste(opt$infile, " N = ", nrow(tbl), " event = ", length(which(tbl$status == 1)), sep =
      "")
  } else {
    as.character(row[1])
  }

  grid.arrange(
    pl_km,
    pl_roc,
    ncol = 2,
    nrow = 1,
    top = textGrob(title,
                   gp = gpar(fontsize = 14))
  )
  dev.off()
}

for (i in rownames(formula)) {
  do_plot(tbl, formula[i,])
}