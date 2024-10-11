rm(list = ls())
setwd("~/fat/amp_bsseq/")
library(reshape2)
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(survival)
library(survminer)
library(timeROC)
library(scales)
library(grid)
library(gridExtra)
library(extrafont)
library(knitr)
library(tibble)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)

plot_km <- function(tbl, fit, diff) {
  kmp <- pchisq(diff$chisq, length(diff$n) - 1, lower.tail = F)
  kmp <- signif(kmp, digits = 2)

  hr <- (diff$obs[1] / diff$exp[1]) / (diff$obs[2] / diff$exp[2])
  cilow <- exp(log(hr) - qnorm(0.975) * sqrt(1 / diff$exp[1] + 1 / diff$exp[2]))
  cihigh <- exp(log(hr) + qnorm(0.975) * sqrt(1 / diff$exp[1] + 1 / diff$exp[2]))
  hr <- signif(hr, digits = 2)
  cilow <- signif(cilow, digits = 2)
  cihigh <- signif(cihigh, digits = 2)

  ggsurv <- ggsurvplot(
    fit,
    data = tbl,
    legend.labs = c(paste("high:", sum(
      tbl$group == "Group 1"
    )),
                    paste("low:", sum(
                      tbl$group == "Group 2"
                    ))),
    linetype = c("solid", "solid"),
    size = 0.8,
    # line width
    palette = c("#b53c32", "#2b628e"),
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
        label = paste0("P = ", kmp, "\nHR = ", hr, "\n95% CI: ", cilow, "-", cihigh),
        x = 0,
        y = 0
      ),
      hjust = 0,
      vjust = 0,
      size = 3
    )

  breaks <- seq(0, 1, 0.2)
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
  if (kmp < 0.05) {
    plot <-
      plot + theme(panel.background = element_rect(fill = "#f7f2e6"))
  }
  return(plot)
}

plot_km_ss <- function(tbl, fit, diff) {
  kmp <- pchisq(diff$chisq, length(diff$n) - 1, lower.tail = F)
  kmp <- signif(kmp, digits = 2)

  hr <- (diff$obs[1] / diff$exp[1]) / (diff$obs[2] / diff$exp[2])
  cilow <- exp(log(hr) - qnorm(0.975) * sqrt(1 / diff$exp[1] + 1 / diff$exp[2]))
  cihigh <- exp(log(hr) + qnorm(0.975) * sqrt(1 / diff$exp[1] + 1 / diff$exp[2]))
  hr <- signif(hr, digits = 2)
  cilow <- signif(cilow, digits = 2)
  cihigh <- signif(cihigh, digits = 2)

  ggsurv <- ggsurvplot(
    fit,
    data = tbl,
    legend.labs = c(paste("steep:", sum(
      tbl$group == "Group 1"
    )),
                    paste("shallow:", sum(
                      tbl$group == "Group 2"
                    ))),
    linetype = c("solid", "solid"),
    size = 0.8,
    # line width
    palette = c("#b53c32", "#2b628e"),
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
        label = paste0("P = ", kmp, "\nHR = ", hr, "\n95% CI: ", cilow, "-", cihigh),
        x = 0,
        y = 0
      ),
      hjust = 0,
      vjust = 0,
      size = 3
    )

  breaks <- seq(0, 1, 0.2)
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
  if (kmp < 0.05) {
    plot <-
      plot + theme(panel.background = element_rect(fill = "#f7f2e6"))
  }
  return(plot)
}

data <- read.table(file = "sites.tsv", sep = "\t", header = T)
data_long <-
  melt(
    data,
    id.vars = c("gene", "chr", "site", "trend"),
    variable.name = "sample",
    value.name = "beta"
  )

data_long <-
  separate(
    data_long,
    col = "sample",
    into = c("sampleID", "type"),
    sep = "_"
  )

data_long <- data_long %>%
  mutate(group_id = paste(gene, chr, site, sep = "_"))

plot_list <- list()

for (group in unique(data_long$group_id)) {
  subset_data <- data_long %>% filter(group_id == group)
  subset_wide <- subset_data %>%
    pivot_wider(names_from = type, values_from = beta)
  subset_wide <- subset_wide %>% filter(!is.na(C) & !is.na(N))
  n_higher <- sum(subset_wide$N > subset_wide$C, na.rm = TRUE)
  c_higher <- sum(subset_wide$C > subset_wide$N, na.rm = TRUE)

  avg_values <- subset_data %>%
    group_by(type) %>%
    summarize(mean_beta = mean(beta, na.rm = TRUE))

  p <-
    ggplot(subset_data,
           aes(
             x = type,
             y = beta,
             group = sampleID,
             color = sampleID
           )) +
      geom_line() +
      geom_point() +
      geom_text(
        data = avg_values,
        aes(
          x = type,
          y = mean_beta,
          label = round(mean_beta, 2)
        ),
        inherit.aes = FALSE,
        vjust = -1,
        color = "black"
      ) +
      labs(
        title = paste(
          unique(subset_data$gene),
          unique(subset_data$chr),
          unique(subset_data$site),
          unique(subset_data$trend)
        ),
        subtitle = paste("N > C:", n_higher, "| C > N:", c_higher),
        x = "Type",
        y = "Beta Value"
      ) +
      theme_minimal()
  plot_list[[group]] <- p
}

pdf("output_plots.pdf", width = 7, height = 5)
for (group in names(plot_list)) {
  print(plot_list[[group]])
}
dev.off()

data_long <- data %>%
  mutate(group_id = paste(gene, chr, site, trend, sep = "_")) %>%
  mutate(group_id = str_replace(group_id, "_$", ""))
data_long <- data_long[, -1:-4]
data_long <- data_long %>%
  pivot_longer(
    cols = starts_with("FXH"),
    names_to = c("sampleID", "type"),
    names_sep = "_",
    values_to = "beta"
  )

data_wide <- data_long %>%
  pivot_wider(names_from = type,
              values_from = beta)
data_wide1 <- data_wide %>%
  mutate(value = (C - N) / C)
data_wide1 <- na.omit(data_wide1)
data_wide1 <- data_wide1[, -3:-4]
data_wide1 <- data_wide1 %>%
  pivot_wider(names_from = sampleID,
              values_from = value)
data_wide1 <- as.data.frame(data_wide1)
row.names(data_wide1) <- data_wide1$group_id
data_wide1 <- data_wide1[, -1]
data_wide1 <- as.data.frame(t(data_wide1))

sampleinfo <-
  read.table(
    "../rawdata/score5.tsv",
    header = T,
    sep = "\t",
    row.names = 1
  )

sampleinfo$`#sample` <- rownames(sampleinfo)
data_wide1$`#sample` <- rownames(data_wide1)
df_combined <-
  merge(sampleinfo[, c(2, 1, ncol(sampleinfo))], data_wide1, by = "#sample", sort = FALSE)
write.table(
  df_combined,
  row.names = F,
  file = "data.tsv",
  sep = "\t",
  quote = F
)

threshold <- 1825
tbl <- read_tsv("data.tsv", col_type = cols())
tbl$status <- ifelse(tbl$time >= threshold, 0, tbl$status)
for (i in 4:63) {
  markers <- colnames(tbl)[i]
  tbl1 <- na.omit(tbl[, c(1, 2, 3, i)])
  tbl1 <- tbl1 %>%
    mutate(
      group = case_when(
        tbl1[[4]] > 0.05 ~ "Group 1",
        tbl1[[4]] <= 0.05 & tbl1[[4]] > -0.05 ~ "Group 2",
        tbl1[[4]] <= -0.05 ~ "Group 1",
        TRUE ~ "Unclassified"
      )
    )
  figfile <- str_interp("./${markers}.pdf")
  write(figfile, stderr())
  table(tbl1$group)
  if (all(table(tbl1$group) > 2) &&
    length(table(tbl1$group)) > 1) {
    diff <- tryCatch({
      survdiff(Surv(time, status) ~ group,
               data = tbl1,
               rho = 0)
    },
      error = function(cond) {
        message(cond, "\n")
        return(NA)
      })

    fit <- tryCatch({
      survfit(Surv(time, status) ~ group, data = tbl1)
    },
      error = function(cond) {
        message(cond, "\n")
        return(NA)
      })
    pl_km <- plot_km(tbl1, fit, diff)
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
  }
  pdf(
    file = figfile,
    family = "Helvetica",
    width = 4,
    height = 4,
    useDingbats = FALSE
  )
  title <-
    paste0(markers, " N = ", nrow(tbl1), " event = ", length(which(tbl1$status == 1)))
  grid.arrange(
    pl_km,
    ncol = 1,
    nrow = 1,
    top = textGrob(title,
                   gp = gpar(fontsize = 10))
  )
  dev.off()
}

dtfr <- read.csv("score.tsv",
                 header = T,
                 sep = "\t",
                 row.names = 1)

HR1 <- function(DF, LINE) {
  dtfr_omit <- na.omit(DF)
  colnames(dtfr_omit) <- c("Status", "Time", "Score")
  dtfr_omit$ScoreAboveLINE <-
    factor(
      ifelse(dtfr_omit$Score > LINE, 0, 1),
      levels = c(0, 1),
      labels = c("Low", "High")
    )
  if (all(table(dtfr_omit$ScoreAboveLINE) > 5) &&
    length(table(dtfr_omit$ScoreAboveLINE)) > 1) {
    hr <-
      coxph(Surv(Time, Status) ~ ScoreAboveLINE, data = dtfr_omit) %>% summary()
    p <-
      survdiff(Surv(Time, Status) ~ ScoreAboveLINE, data = dtfr_omit)
    return(list(
      p_value = p$pvalue,
      hr_coefficient = hr$coefficients[2]
    ))
  } else {
    return(list(p_value = NA,
                hr_coefficient = NA))
  }
}

myHR <- list()
myP <- list()
pdf("km_plots.pdf", width = 7, height = 5)
for (i in 3:62) {
  dtfr0 <- dtfr[, c(1, 2, i)]
  title <- colnames(dtfr)[i]
  colnames(dtfr0) <- c("Status", "Time", "Score")
  LINE1 <- 0.05
  result <- HR1(dtfr0, LINE1)
  myHR$ratio[i - 2] <- result$hr_coefficient
  myP$ratio[i - 2] <- result$p_value
  dtfr_omit <- na.omit(dtfr0)
  colnames(dtfr_omit) <- c("Status", "Time", "Score")
  dtfr_omit$ScoreAboveLINE <-
    factor(
      ifelse(dtfr_omit$Score > LINE1, 0, 1),
      levels = c(0, 1),
      labels = c("Low", "High")
    )
  p <-
    survfit2(Surv(Time, Status) ~ ScoreAboveLINE, data = dtfr_omit) %>%
      ggsurvfit() +
      labs(title = paste(title, "OS"),
           x = "Days",
           y = "Overall survival probability") +
      add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
      add_pvalue(location = "annotation") +
      scale_ggsurvfit()
  print(ggsurvfit_build(p))
}
dev.off()
