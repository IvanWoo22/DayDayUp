#!/usr/bin/env Rscript
library(pROC, quietly = T)
library(data.table, quietly = T)
library(readr, quietly = T)

Args <- commandArgs(T)
target <- Args[1]
data_path1 <- Args[2]
data_path2 <- Args[3]
output_path <- Args[4]

do_logistic_ML <-
  function(df_train,
           df_test,
           ProbeID) {
    formula_string <- paste("judge ~ ", ProbeID, sep = "")
    df_train <- na.omit(df_train)
    df_test <- na.omit(df_test)
    median_diff <-
      abs(median(df_train[df_train$judge == T, ProbeID]) - median(df_train[df_train$judge ==
                                                                             F, ProbeID]))
    mean_diff <-
      abs(mean(df_train[df_train$judge == T, ProbeID]) - mean(df_train[df_train$judge ==
                                                                         F, ProbeID]))
    res_logistic <-
      glm(as.formula(formula_string), df_train, family = "binomial")
    
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
    
    testauc <- format(round(
      roc(
        df_test$judge,
        predict(res_logistic, df_test, type = "link"),
        levels = c(F, T),
        direction = "<"
      )$auc,
      5
    ), nsmall = 5)
    
    if (nrow(summary(res_logistic)$coefficients) > 1) {
      return(as.data.frame(
        cbind(
          ProbeID,
          median_diff,
          mean_diff,
          format(round(
            summary(res_logistic)$coefficients[2, 1], 5
          ), nsmall = 5),
          format(round(
            summary(res_logistic)$coefficients[2, 4], 5
          ), nsmall = 5),
          format(round(((
            summary(res_logistic)$null.deviance / -2
          ) - (summary(res_logistic)$deviance / -2)) / (summary(res_logistic)$null.deviance / -2), 5
          ), nsmall = 5),
          format(round(1 - pchisq(2 * (
            (summary(res_logistic)$deviance / -2) - (summary(res_logistic)$null.deviance / -2)
          ), df = 1), 10), nsmall = 10),
          rocauc,
          testauc
        ),
        stringAsFactor = F
      ))
    } else{
      return(as.data.frame(
        cbind(
          ProbeID,
          median_diff,
          mean_diff,
          "NA",
          "NA",
          "NA",
          "NA",
          rocauc,
          testauc
        ),
        stringAsFactor = F
      ))
    }
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
train_df[train_df$status == target,]$judge <- T
train_df$judge <- as.factor(train_df$judge)
test_df$judge <- F
test_df[test_df$status == target,]$judge <- T
test_df$judge <- as.factor(test_df$judge)

write_delim(
  data.frame(t(
    c(
      "#marker",
      "median_diff",
      "mean_diff",
      "coef",
      "coef_p",
      "R_2",
      "reg_p",
      "rocauc",
      "testauc"
    )
  )),
  output_path,
  append = TRUE,
  col_names = F,
  delim = "\t"
)

for (ProbeID in colnames(train_df)[3:ncol(train_df) - 1]) {
  output <-
    do_logistic_ML(train_df[, which(names(train_df) %in% c(ProbeID, "judge"))],
                   test_df[, which(names(test_df) %in% c(ProbeID, "judge"))],
                   ProbeID)
  write_delim(
    output,
    output_path,
    append = TRUE,
    col_names = F,
    delim = "\t"
  )
}