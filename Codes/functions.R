# CV lasso and mod LASSO


# Function for summary statistics lasso
summary.lasso <- function(lasso.mod, x, y, ind) {
  
  x.train <- x[rownames(x) %in% ind, ]
  x.test <- x[!rownames(x) %in% ind, ]
  
  # Prediction
  pred.train <- try(as.factor(predict(lasso.mod, newx=x.train, type="class",
                                      s="lambda.min", alpha=1)), silent=TRUE)
  pred.test <- try(as.factor(predict(lasso.mod, newx=x.test, type="class",
                                     s="lambda.min", alpha=1)), silent = TRUE)
  
  # Confusion matrix
  conf.train <- caret::confusionMatrix(pred.train, as.factor(y[names(y) %in% ind]))
  conf.test <- caret::confusionMatrix(pred.test, as.factor(y[!names(y) %in% ind]))
  
  # Summary statistics
  df.eval <- data.frame(Type=c("Train", "Test"),
                        Kappa=c(conf.train$overall[["Kappa"]],
                                conf.test$overall[["Kappa"]]
                        ),
                        `Balanced Acc`=c(conf.train$byClass[["Balanced Accuracy"]],
                                         conf.test$byClass[["Balanced Accuracy"]]
                        ),
                        Sensitivity=c(conf.train$byClass[["Sensitivity"]],
                                      conf.test$byClass[["Sensitivity"]]
                        ),
                        Specificity=c(conf.train$byClass[["Specificity"]],
                                      conf.test$byClass[["Specificity"]]
                        ),
                        PPV=c(conf.train$byClass[["Pos Pred Value"]],
                              conf.test$byClass[["Pos Pred Value"]]
                        ),
                        NPV=c(conf.train$byClass[["Neg Pred Value"]],
                              conf.test$byClass[["Neg Pred Value"]]
                        )
  )
  
  
  # Selected variable
  selected.var <- as.data.frame.matrix(coef(lasso.mod))
  names(selected.var) <- 'coef'
  selected.var$var <- rownames(selected.var)
  rownames(selected.var) <- NULL
  selected.var <- selected.var[selected.var$coef!=0 & selected.var$var!='(Intercept)', ]
  selected.var <- selected.var[, c(2, 1)]
  
  
  # Output
  out <- list("evaluation"=df.eval, 
              "pred.test"=pred.test, 
              'selected.var'=selected.var
  )
}


summary.glm <- function(glm.mod, x, y, ind) {
  
  if(!is.data.frame(x)) {
    new.x.train <- as.data.frame(x)
  }
  
  new.x.train <- x[rownames(x) %in% ind, ]
  new.x.test <- x[!rownames(x) %in% ind, ]
  y.train <- y[names(y) %in% ind]
  y.test <- y[!names(y) %in% ind]
  
  # Prediction
  pred.train <- predict(glm.mod, newx=new.x.train, type="response")
  pred.train <- as.factor(ifelse(pred.train > 0.5, 'CRC', 'Adenoma'))
  
  pred.test <- predict(glm.mod, newdata=as.data.frame(new.x.test), type="response")
  pred.test <- as.factor(ifelse(pred.test > 0.5, 'CRC', 'Adenoma'))
  
  
  # Confusion matrix
  conf.train <- caret::confusionMatrix(pred.train, as.factor(y.train))
  conf.test <- caret::confusionMatrix(pred.test, as.factor(y.test))
  
  # Summary statistics
  df.eval <- data.frame(Type= c("Train", "Test"),
                        Kappa=c(conf.train$overall[["Kappa"]],
                                conf.test$overall[["Kappa"]]
                        ),
                        `Balanced Acc`=c(conf.train$byClass[["Balanced Accuracy"]],
                                         conf.test$byClass[["Balanced Accuracy"]]
                        ),
                        Sensitivity=c(conf.train$byClass[["Sensitivity"]],
                                      conf.test$byClass[["Sensitivity"]]
                        ),
                        Specificity=c(conf.train$byClass[["Specificity"]],
                                      conf.test$byClass[["Specificity"]]
                        ),
                        PPV=c(conf.train$byClass[["Pos Pred Value"]],
                              conf.test$byClass[["Pos Pred Value"]]
                        ),
                        NPV=c(conf.train$byClass[["Neg Pred Value"]],
                              conf.test$byClass[["Neg Pred Value"]]
                        )
  )
  
  
  # Output
  out <- list("evaluation"=df.eval, 
              "pred.test"=pred.test
  )
}
