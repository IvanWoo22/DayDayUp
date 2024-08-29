rm(list = ls())
library(ggplot2)
library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(gridExtra)
library(grid)
setwd("~/fat/rawdata/")

dtfr <- read.csv("score4.tsv",
                 header = T,
                 sep = "\t",
                 row.names = 1)

dtfr_omit <- na.omit(dtfr)
dtfr$Status <-
  factor(dtfr$Status,
         levels = c("1", "0"),
         labels = c("dead", "censored"))
dtfr$Med <-
  factor(dtfr$Med,
         levels = c("1", "0"),
         labels = c("治疗", "无治疗"))
ggplot(dtfr, aes(
  x = reorder(Sample, Score),
  y = Score,
  fill = Status
)) +
  geom_bar(stat = "identity") +
  labs(title = "Sample Scores", x = "Sample", y = "Risk Score") +
  coord_flip() +
  theme_minimal()

ggplot(dtfr_omit, aes(
  x = reorder(Sample, Score),
  y = Score,
  fill = Status
)) +
  geom_bar(stat = "identity") +
  labs(title = "Sample Scores", x = "Sample", y = "Risk Score") +
  coord_flip() +
  theme_minimal()

Surv(dtfr$Time, dtfr$Status)[1:10]
head(lung[, c("time", "status", "sex")])
survfit2(Surv(Time, Status) ~ 1, data = dtfr) %>%
  ggsurvfit() +
  labs(title = "Total OS",
       x = "Days",
       y = "Overall survival probability")

survfit2(Surv(Time, Status) ~ Sex, data = dtfr) %>%
  ggsurvfit() +
  labs(title = "Sex OS",
       x = "Days",
       y = "Overall survival probability")
survdiff(Surv(Time, Status) ~ Sex, data = dtfr)

dtfr$AgeAboveMedian <- ifelse(dtfr$Age > median(dtfr$Age), 1, 0)
dtfr$AgeAboveMean <- ifelse(dtfr$Age > mean(dtfr$Age), 1, 0)
survfit2(Surv(Time, Status) ~ AgeAboveMedian, data = dtfr) %>%
  ggsurvfit() +
  labs(title = "AgeAboveMedian OS",
       x = "Days",
       y = "Overall survival probability")
survdiff(Surv(Time, Status) ~ AgeAboveMedian, data = dtfr)

survfit2(Surv(Time, Status) ~ AgeAboveMean, data = dtfr) %>%
  ggsurvfit() +
  labs(title = "AgeAboveMean OS",
       x = "Days",
       y = "Overall survival probability")
survdiff(Surv(Time, Status) ~ AgeAboveMean, data = dtfr)

dtfr_omit_high <-
  dtfr_omit[which(dtfr_omit$Score > median(dtfr$Score, na.rm = T)),]
survfit2(Surv(Time, Status) ~ Med, data = dtfr_omit_high) %>%
  ggsurvfit() +
  labs(title = "Med OS in high risk group",
       x = "Days",
       y = "Overall survival probability")
survdiff(Surv(Time, Status) ~ Med, data = dtfr_omit_high)

dtfr_omit$MedF <-
  factor(dtfr_omit$Med,
         levels = c("1", "0"),
         labels = c("W/ Med", "W/O Med"))
for (i in 8:429) {
  dtfr_omit$ScoreAbove <-
    ifelse(as.array(dtfr_omit[, i]) > median(as.array(dtfr[, i]), na.rm = T),
           "High",
           "Low")
  dtfr_omit1 <- dtfr_omit[dtfr_omit$Med == 1,]
  dtfr_omit2 <- dtfr_omit[dtfr_omit$Med == 0,]
  dtfr_omit3 <-
    dtfr_omit[which(as.array(dtfr_omit[, i]) > median(as.array(dtfr[, i]), na.rm = T)),]
  dtfr_omit4 <-
    dtfr_omit[which(as.array(dtfr_omit[, i]) <= median(as.array(dtfr[, i]), na.rm = T)),]

  if (all(table(dtfr_omit1$ScoreAbove) > 5) &
    all(table(dtfr_omit3$Med) > 5) &
    all(table(dtfr_omit4$Med) > 5)) {
    p1 <-
      survfit2(Surv(Time, Status) ~ ScoreAbove, data = dtfr_omit1) %>%
        ggsurvfit() +
        labs(title = "W/ Med OS",
             x = "Days",
             y = "Overall survival probability") +
        add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
        add_pvalue(location = "annotation") +
        scale_ggsurvfit()
    p2 <-
      survfit2(Surv(Time, Status) ~ ScoreAbove, data = dtfr_omit2) %>%
        ggsurvfit() +
        labs(title = "W/O Med OS",
             x = "Days",
             y = "Overall survival probability") +
        add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
        add_pvalue(location = "annotation") +
        scale_ggsurvfit()
    p3 <-
      survfit2(Surv(Time, Status) ~ MedF, data = dtfr_omit3) %>%
        ggsurvfit() +
        labs(title = "High Risk OS",
             x = "Days",
             y = "Overall survival probability") +
        add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
        add_pvalue(location = "annotation") +
        scale_ggsurvfit()
    p4 <-
      survfit2(Surv(Time, Status) ~ MedF, data = dtfr_omit4) %>%
        ggsurvfit() +
        labs(title = "Low Risk OS",
             x = "Days",
             y = "Overall survival probability") +
        add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
        add_pvalue(location = "annotation") +
        scale_ggsurvfit()
    outtext <- paste("Score", i - 7, ".pdf", sep = "")
    titletext <- paste("Score", i - 7, sep = "")
    pdf(outtext, width = 8.7, height = 8.7)
    grid.arrange(p1,
                 p2,
                 p3,
                 p4,
                 ncol = 2,
                 top = textGrob(titletext, gp = gpar(
                   fontsize = 14, font = 3
                 )))
    dev.off()
  }
}

for (i in 10:432) {
  dtfr1 <- dtfr[, c(1, 2, i)]
  dtfr_omit <- na.omit(dtfr1)
  colnames(dtfr_omit) <- c("Status", "Time", "Score")
  LINE <-
    dtfr_omit$Score[order(dtfr_omit$Score, decreasing = TRUE)[sum(dtfr_omit$Status, na.rm = T)]]

  dtfr_omit$ScoreAboveMedian <-
    ifelse(dtfr_omit$Score > median(dtfr_omit$Score), "High", "Low")
  p1 <-
    survfit2(Surv(Time, Status) ~ ScoreAboveMedian, data = dtfr_omit) %>%
      ggsurvfit() +
      labs(title = "ScoreAboveMedian OS",
           x = "Days",
           y = "Overall survival probability") +
      add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
      add_pvalue(location = "annotation") +
      scale_ggsurvfit()

  dtfr_omit$ScoreAboveMean <-
    ifelse(dtfr_omit$Score > mean(dtfr_omit$Score), "High", "Low")
  p2 <-
    survfit2(Surv(Time, Status) ~ ScoreAboveMean, data = dtfr_omit) %>%
      ggsurvfit() +
      labs(title = "ScoreAboveMean OS",
           x = "Days",
           y = "Overall survival probability") +
      add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
      add_pvalue(location = "annotation") +
      scale_ggsurvfit()

  dtfr_omit$ScoreAboveLINE <-
    ifelse(dtfr_omit$Score > LINE, "High", "Low")
  p3 <-
    survfit2(Surv(Time, Status) ~ ScoreAboveLINE, data = dtfr_omit) %>%
      ggsurvfit() +
      labs(title = "ScoreAboveLINE OS",
           x = "Days",
           y = "Overall survival probability") +
      add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
      add_pvalue(location = "annotation") +
      scale_ggsurvfit()
  outtext <- paste("TScore", i - 9, ".pdf", sep = "")
  titletext <- paste("Score", i - 9, sep = "")
  pdf(outtext, width = 11.7, height = 4.7)
  grid.arrange(p1,
               p2,
               p3,
               ncol = 3,
               top = textGrob(titletext, gp = gpar(fontsize = 14, font = 3)))
  dev.off()
}

i = 333
dtfr0 <- dtfr[, c(1, 2, i)]
colnames(dtfr0) <- c("Status", "Time", "Score")
LINE <-
  dtfr0$Score[order(dtfr0$Score, decreasing = TRUE)[sum(dtfr0$Status, na.rm = T)]]

dtfr1 <- dtfr[which(dtfr$Sex == "Male"), c(1, 2, i)]
dtfr_omit <- na.omit(dtfr1)
colnames(dtfr_omit) <- c("Status", "Time", "Score")

dtfr_omit$ScoreAboveLINE <-
  ifelse(dtfr_omit$Score > LINE, "High", "Low")
survfit2(Surv(Time, Status) ~ ScoreAboveLINE, data = dtfr_omit) %>%
  ggsurvfit() +
  labs(title = "ScoreAboveLINE OS",
       x = "Days",
       y = "Overall survival probability") +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  add_pvalue(location = "annotation") +
  scale_ggsurvfit()

dtfr2 <- dtfr[which(dtfr$Sex == "Female"), c(1, 2, i)]
dtfr_omit <- na.omit(dtfr2)
colnames(dtfr_omit) <- c("Status", "Time", "Score")

dtfr_omit$ScoreAboveLINE <-
  factor(
    ifelse(dtfr_omit$Score > LINE, 1, 0),
    levels = c(0, 1),
    labels = c("Low", "High")
  )

survfit2(Surv(Time, Status) ~ ScoreAboveLINE, data = dtfr_omit) %>%
  ggsurvfit() +
  labs(title = "ScoreAboveLINE OS",
       x = "Days",
       y = "Overall survival probability") +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  add_pvalue(location = "annotation") +
  scale_ggsurvfit()
t <- survdiff(Surv(Time, Status) ~ ScoreAboveLINE, data = dtfr_omit)
t$pvalue

coxph(Surv(Time, Status) ~ ScoreAboveLINE, data = dtfr_omit) %>%
  summary()
survdiff(Surv(Time, Status) ~ ScoreAboveLINE, data = dtfr_omit) %>%
  summary()
