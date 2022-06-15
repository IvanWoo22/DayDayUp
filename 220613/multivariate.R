#!/usr/bin/env Rscript
library(pROC)
library(data.table)
library(readr)

Args <- commandArgs(T)
target <- Args[1]
data_path1 <- Args[2]
data_path2 <- Args[3]
feature_path <- Args[4]
output_path <- Args[5]

do_logistic_ML <- function(df_train, df_test, df_feature) {
  formula_string <- paste("judge ~ ", df_feature, sep = "")
  res_logistic <-
    glm(as.formula(formula_string), df_train, family = "binomial")
  
  #------------------------------------------------------#
  # ROC
  #------------------------------------------------------#
  if (length(df_train$judge) == length(res_logistic$fitted.values)) {
    rocauc <- format(round(
      roc(
        df_train$judge,
        res_logistic$fitted.values,
        levels = c(F, T),
        direction = "<"
      )$auc,
      5
    ), nsmall = 5)
  } else {
    rocauc <- NA
  }
  
  #------------------------------------------------------#
  # ML test AUC value
  #------------------------------------------------------#
  testauc <- format(round(
    roc(
      df_test$judge,
      predict(res_logistic, df_test, type = "link"),
      levels = c(F, T),
      direction = "<"
    )$auc,
    5
  ), nsmall = 5)
  
  # Output
  return(as.data.frame(
    cbind(
      df_feature,
      as.numeric(format(
        round(summary(res_logistic)$coefficients[2, 1], 5), nsmall = 5
      )),
      as.numeric(format(
        round(summary(res_logistic)$coefficients[2, 4], 5), nsmall = 5
      )),
      format(round(((summary(res_logistic)$null.deviance / -2) - (summary(res_logistic)$deviance / -2)
      ) / (
        summary(res_logistic)$null.deviance / -2
      ), 5), nsmall = 5),
      format(round(1 - pchisq(2 * (
        (summary(res_logistic)$deviance / -2) - (summary(res_logistic)$null.deviance / -2)
      ), df = 1), 10), nsmall = 10),
      rocauc,
      testauc
    ),
    stringAsFactor = F
  ))
}

train_df <-
  read.csv(data_path1,
           header = T,
           row.names = 1,
           sep = "\t")
test_df <-
  read.csv(data_path2,
           header = T,
           row.names = 1,
           sep = "\t")
train_df$judge <- F
train_df[train_df$status == target, ]$judge <- T
train_df$judge <- as.factor(train_df$judge)
test_df$judge <- F
test_df[test_df$status == target, ]$judge <- T
test_df$judge <- as.factor(test_df$judge)

feature <-
  as.data.frame(fread(feature_path, header = T), stringAsFactor = F)

write_delim(
  data.frame(t(
    c("ProbeID", "coef", "coef_p", "R_2", "reg_p", "rocauc", "testauc")
  )),
  output_path,
  append = TRUE,
  col_names = F,
  delim = "\t"
)

for (i in 1:nrow(feature)) {
  output <-
    do_logistic_ML(train_df, test_df, feature[i, 1])
  write_delim(
    output,
    output_path,
    append = TRUE,
    col_names = F,
    delim = "\t"
  )
}