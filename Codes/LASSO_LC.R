##########################################################################################
## PERMUTATION HIGH DIMENSIONAL DATA                                                    ##
## Data         : lung cancer data                                                      ##
## Prepared by  : Dea Putri (dputri@its.jnj.com)                                        ##
##########################################################################################
load(file = "Data/ALLDATA.RData")

clin <- c("Group", "Age", "Packyears","BMI", 
          "COPD", "Cardiac", "Smoking", "Coagulation")
meta <- names(ALL.train)[grepl("var", names(ALL.train), ignore.case = TRUE)]
var <- c(clin, meta)

ALL.train[, "Smoking"] <- factor(ALL.train[, "Smoking"],
                                 levels = c("Never", "Active", "Stopped"))
table(ALL.train[, "Smoking"])
ALL.test[, "Smoking"] <- factor(ALL.test[, "Smoking"],
                                levels = c("Never", "Active", "Stopped"))
table(ALL.test[, "Smoking"])


# Create the matrix for lasso input
y <- ALL.train[, "Group"]
X <- ALL.train[var]
X <- subset(X, select = -c(Group))
X <- data.matrix(X)

penalty.fctr <- rep(0, ncol(X))
names(penalty.fctr) <- colnames(X)
penalty.fctr[meta] <- 1


# Create response for test data
X.test <- ALL.test[var]
X.test <- subset(X.test, select = -c(Group))
X.test <- data.matrix(X.test)

y.test <- ALL.test[, "Group"]





#### LASSO 
fit.lasso <- function(X, y, alpha = 1, partition = 0.7, penalty = NULL, new.x, new.y, n.iter = 1000)  {
  require(glmnet)
  require(caret)
  
  out <- vector("list", n.iter)
  
  for(i.k in 1:n.iter) {
    options(warn = -1)
    
    # Create data partition
    ind <- createDataPartition(y, times = 1, p = partition, list = FALSE)
    
    if(is.null(penalty)) {
      cv.lasso <- try(cv.glmnet(x = X[ind, ], y = y[ind], alpha = alpha,
                                standardize = FALSE, family = "binomial",
                                nfolds = 3, type.measure = "class"),
                      silent = TRUE)
    } else {
      cv.lasso <- try(cv.glmnet(x = X[ind, ], y = y[ind], alpha = alpha,
                                standardize = FALSE, family = "binomial",
                                penalty.factor = penalty, 
                                nfolds = 3, type.measure = "class"),
                      silent = TRUE)
    }
    
    pred.train <- try(as.factor(predict(cv.lasso, newx = X[ind, ], type = "class")), silent = TRUE)
    pred.test <- try(as.factor(predict(cv.lasso, newx = X[-ind, ], type = "class")), silent = TRUE)
    pred.val <- try(as.factor(predict(cv.lasso, newx = new.x, type = "class")), silent = TRUE)
    
    conf.train <- caret::confusionMatrix(pred.train, as.factor(y[ind]))
    conf.test <- caret::confusionMatrix(pred.test, as.factor(y[-ind]))
    conf.val <- caret::confusionMatrix(pred.val, as.factor(new.y))
    
    df.eval <- data.frame(Type = c("Train", "Test", "Validation"),
                          Kappa = c(conf.train$overall[["Kappa"]],
                                    conf.test$overall[["Kappa"]],
                                    conf.val$overall[["Kappa"]]),
                          `Balanced Acc` =  c(conf.train$byClass[["Balanced Accuracy"]],
                                              conf.test$byClass[["Balanced Accuracy"]],
                                              conf.val$byClass[["Balanced Accuracy"]])
    )
    
    train.ind <- ind
    
    if(n.iter == 1) {
      out <- list("evaluation" = df.eval, 
                  "train.ind" = train.ind, 
                  "pred.train" = pred.train, 
                  "pred.test" = pred.test, 
                  "pred.val" = pred.val)
    } else {
      out[[i.k]] <- list("evaluation" = df.eval, 
                         "train.ind" = train.ind, 
                         "pred.train" = pred.train, 
                         "pred.test" = pred.test, 
                         "pred.val" = pred.val)
      names(out)[[i.k]] <- i.k
      print(i.k)
    }
  }
  return(out)
}


test <- fit.lasso(X = X, y = y, penalty = penalty.fctr, new.x = X.test, new.y = y.test, n.iter = 1)

# Not-standardized
set.seed(1)
t2 = Sys.time()
lc.lasso <- fit.lasso(X = X, y = y, penalty = penalty.fctr, new.x = X.test, new.y = y.test)
t <- Sys.time()
print((t - t2))    # 53.29 mins
save(lc.lasso, file = "Data/lc_lasso.Rda")





#### LASSO: PERMUTATION
clin <- c("Age", "Packyears","BMI", 
          "COPD", "Cardiac", "Smoking", "Coagulation")

n.iter <- 1000
out <- vector("list", n.iter)


t2 = Sys.time()
for(i.k in 1:n.iter) {
  set.seed(i.k)
  smpl <- sample(nrow(X))
  X.train.new <- data.frame(X[, clin], X[smpl, meta])
  X.train.new <- data.matrix(X.train.new)
  
  out[[i.k]] <- try(fit.lasso(X = X.train.new, y = y, penalty = penalty.fctr,
                              new.x = X.test, new.y = y.test, n.iter = 1))
  names(out)[[i.k]] <- i.k
  print(i.k)
}
t <- Sys.time()
print((t - t2))    # 2.05 hours
lc.lasso.perm <- out
save(lc.lasso.perm, file = "Data/lc_lasso_perm.Rda")


