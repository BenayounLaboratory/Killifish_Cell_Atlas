options(stringsAsFactors = FALSE)
library(pROC)
library(caret)
library(randomForest)    # random forest

####################  Aaron functions
getTwoClassBalancedAccuracy <- function (data, lev = NULL, model = NULL) {
  require(pROC)
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))){
    stop("levels of observed and predicted data do not match")
  }
  rocObject <- try(pROC:::roc(data$obs, data[, lev[1]]), silent = TRUE)
  if (class(rocObject)[1] == "try-error") {
    return(NA)
  }else{
    out <- ( sensitivity(data[, "pred"], data[, "obs"], lev[1]) + specificity(data[, "pred"], data[, "obs"], lev[2]) ) /2
    names(out) <- "balancedAcc"
    return(out)
  }
}




