#!/usr/bin/env Rscript
library(pROC, quietly = T)
library(data.table, quietly = T)
library(readr, quietly = T)

Args <- commandArgs(T)
target <- Args[1]
data_path <- Args[2]
feature_path <- Args[3]
output_path <- Args[4]

do_logistic_ML <- function(df_test, ProbeID) {
  formula_string <- paste("judge ~ ", ProbeID, sep = "")
  res_logistic <-
    glm(as.formula(formula_string), df_test, family = "binomial")
  
  #------------------------------------------------------#
  # ROC
  #------------------------------------------------------#
  if (length(df_test$judge) == length(res_logistic$fitted.values)) {
    rocauc <- format(round(
      roc(
        df_test$judge,
        res_logistic$fitted.values,
        levels = c(F, T),
        direction = "<"
      )$auc,
      5
    ), nsmall = 5)
  } else {
    rocauc <- NA
  }
  
  return(as.data.frame(
    cbind(
      ProbeID,
      format(round(
        summary(res_logistic)$coefficients[2, 1], 5
      ), nsmall = 5),
      format(round(
        summary(res_logistic)$coefficients[2, 4], 5
      ), nsmall = 5),
      format(round(((summary(res_logistic)$null.deviance / -2) - (summary(res_logistic)$deviance / -2)
      ) / (
        summary(res_logistic)$null.deviance / -2
      ), 5), nsmall = 5),
      format(round(1 - pchisq(2 * (
        (summary(res_logistic)$deviance / -2) - (summary(res_logistic)$null.deviance / -2)
      ), df = 1), 10), nsmall = 10),
      rocauc
    ),
    stringAsFactor = F
  ))
}

test_df <-
  read.csv(data_path,
           header = T,
           sep = "\t")

test_df$judge <- F
test_df[test_df$status == target, ]$judge <- T
test_df$judge <- as.factor(test_df$judge)
feature <-
  as.data.frame(fread(feature_path, header = T), stringAsFactor = F)

write_delim(
  data.frame(t(
    c("ProbeID", "coef", "coef_p", "R_2", "reg_p", "rocauc")
  )),
  output_path,
  append = TRUE,
  col_names = F,
  delim = "\t"
)

for (i in 1:nrow(feature)) {
  output <-
    do_logistic_ML(test_df, feature[i, 1])
  write_delim(
    output,
    output_path,
    append = TRUE,
    col_names = F,
    delim = "\t"
  )
}
