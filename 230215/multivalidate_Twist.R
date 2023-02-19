#!/usr/bin/env Rscript
library(pROC, quietly = T)
library(data.table, quietly = T)
library(readr, quietly = T)

Args <- commandArgs(T)
target <- Args[1]
data_path <- Args[2]
feature_path <- Args[3]
output_path <- Args[4]

test_df <-
  read.csv(data_path,
           header = T,
           sep = "\t")

test_df$judge <- F
test_df[test_df$status == target,]$judge <- T
test_df$judge <- as.factor(test_df$judge)
feature <-
  as.data.frame(fread(feature_path, header = T), stringAsFactor = F)

write_delim(
  data.frame(t(
    c("#marker", "coef", "coef_p", "R_2", "reg_p", "rocauc")
  )),
  output_path,
  append = TRUE,
  col_names = F,
  delim = "\t"
)

for (i in 1:nrow(feature)) {
  df_test <-
    na.omit(test_df[, which(names(test_df) %in% c(strsplit(feature[i, 1], split = "\\+")[[1]], "judge"))])
  if ((nrow(df_test) > 1) &&
      (table(df_test$judge)[1] * table(df_test$judge)[2] > 0)) {
    rocauc <- format(round(
      roc(
        df_test$judge,
        as.data.frame(as.matrix(cbind(1, df_test[, which(names(df_test) %in% c(strsplit(feature[i, 1], split = "\\+")[[1]]))])) %*% as.matrix(as.numeric(
          strsplit(feature[i, 2], split = ",")[[1]]
        )))$V1,
        levels = c(F, T),
        direction = "<"
      )$auc,
      5
    ), nsmall = 5)
    output <- as.data.frame(cbind(feature[i, 1:5],
                                  rocauc),
                            stringAsFactor = F)
  } else {
    output <- as.data.frame(cbind(df_feature,
                                  "NA",
                                  "NA",
                                  "NA",
                                  "NA",
                                  "NA"),
                            stringAsFactor = F)
  }
  write_delim(
    output,
    output_path,
    append = TRUE,
    col_names = F,
    delim = "\t"
  )
}
