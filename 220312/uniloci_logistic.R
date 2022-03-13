#!/usr/bin/env Rscript
library(pROC)
library(data.table)
library(readr)
library(foreach)
library(doParallel)

Args <- commandArgs(T)
data_path <- Args[1]
output_path <- Args[2]

write("==> Config done!", stderr())

do_logistic_ML <- function(df_train, df_test, ProbeID) {
    formula_string <- paste("status ~ ", ProbeID, sep = " ")
    res_logistic <-
        glm(as.formula(formula_string), df_train, family = "binomial")

    if (length(df_train$status) == length(res_logistic$fitted.values)) {
        rocauc <- format(round(
            roc(
                df_train$status,
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
            df_test$status,
            predict(res_logistic, df_test, type = "link"),
            levels = c(F, T),
            direction = "<"
        )$auc,
        5
    ), nsmall = 5)

    return(as.data.frame(
        cbind(
            ProbeID,
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
}

load(data_path)

write("==> Readin over!", stderr())

write_delim(
    data.frame(t(
        c("#ProbeID", "coef", "p", "R_2", "reg_p", "rocauc", "testauc")
    )),
    output_path,
    append = F,
    col_names = F,
    delim = "\t"
)

cl <- makeCluster(8)
registerDoParallel(cl)

write("==> Start parallel!", stderr())

temp <- foreach(
    ProbeID = colnames(train_df)[-ncol(train_df)],
    .packages = c("pROC", "readr"),
    .inorder = F
) %dopar% {
    output <-
        do_logistic_ML(train_df[, which(names(train_df) %in% c(ProbeID, "status"))], test_df[, which(names(test_df) %in% c(ProbeID, "status"))], ProbeID)
    write_delim(
        output,
        output_path,
        append = T,
        col_names = F,
        delim = "\t"
    )
}
quit()