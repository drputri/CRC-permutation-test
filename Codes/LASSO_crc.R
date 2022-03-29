load("Data/group_crcad.Rda")
load("Data/df_crc.Rda")
source("Codes/functions.R")


#####-------------------- CODES --------------------#####
# Create the matrix for lasso input
y <- group
X <- data.matrix(df.crc)

uBiome <- colnames(X)[grepl("SVs", colnames(X), ignore.case = FALSE)]
uArray <- colnames(X)[grepl("_at", colnames(X), ignore.case = FALSE)]

penalty.fctr <- rep(0, ncol(X))
names(penalty.fctr) <- colnames(X)
penalty.fctr[uBiome] <- 1


# LASSO
fit.lasso <- function(X, y, alpha=1, partition=0.6, penalty=NULL, 
                      n.iter=1000) {
  
  require(glmnet)
  require(caret)
  
  out <- vector("list", n.iter)
  
  
  for(i.k in 1:n.iter) {
    options(warn = -1)
    
    
    # Create data partition
    ind <- createDataPartition(y, times = 1, p = partition, list = FALSE)
    
    
    # LOGISTIC REGRESSION
    dt.log <- data.frame(group=y[ind],
                         X[ind, uArray])
    names(dt.log)[2:ncol(dt.log)] <- substring(names(dt.log)[2:ncol(dt.log)], 2)
    mod.glm <- glm(group ~ ., data=dt.log, family=binomial(link='logit'))
    out.glm <- summary.glm(mod.glm, x=as.data.frame(X), y=y, ind=ind) 
    
    
    # LASSO
    if(!is.null(penalty)) {
      cv.lasso <- try(cv.glmnet(x=X[ind, ], y=y[ind], alpha=alpha,
                                standardize=TRUE, family="binomial",
                                penalty.factor=penalty, 
                                nfolds=10, type.measure="class"),
                      silent = TRUE)
      
    } else {
      cv.lasso <- try(cv.glmnet(x=X[ind, ], y=y[ind], alpha=alpha,
                                standardize=TRUE, family="binomial",
                                nfolds=10, parallel=TRUE, type.measure="class"),
                      silent = TRUE)
    }
    
    mod.lasso <- glmnet(x=X[ind, ], y=y[ind], alpha=1, family="binomial",
                        penalty.factor=penalty.fctr, standardize=TRUE,
                        lambda = cv.lasso$lambda.min)
    
    out.lasso <- summary.lasso(mod.lasso, x=X, y=y, ind=ind)
    
    
    # PERMUTATED LASSO
    smpl <- sample(nrow(X))
    X.train.new <- data.frame(X[, uArray], X[smpl, uBiome])
    colnames(X.train.new)[1:14] <- substring(colnames(X.train.new)[1:14], 2)
    X.train.new <- data.matrix(X.train.new)
    
    if(is.null(penalty)) {
      cv.perm <- try(cv.glmnet(x=X.train.new[ind, ], y=y[ind], alpha=alpha,
                               standardize=TRUE, family="binomial",
                               nfolds=10, parallel=TRUE, type.measure="class"),
                     silent = TRUE)
    } else {
      cv.perm <- try(cv.glmnet(x=X.train.new[ind, ], y=y[ind], alpha=alpha,
                               standardize=TRUE, family="binomial",
                               penalty.factor=penalty, 
                               nfolds=10, type.measure="class"),
                     silent = TRUE)
    }
    
    mod.perm <- glmnet(x=X.train.new[ind, ], y=y[ind], alpha=1, family="binomial",
                       penalty.factor=penalty.fctr, standardize=TRUE,
                       lambda = cv.perm$lambda.min)
    
    out.perm <- summary.lasso(mod.perm, x=X.train.new, y=y, ind=ind)
    
    
    # Save the output
    out[[i.k]] <- list(out.glm, out.lasso, out.perm)
    names(out[[i.k]]) <- c("glm", "lasso", "permutation")
    names(out)[[i.k]] <- i.k
    
  }
  
  return(out)
  
}


# Standardized
set.seed(1)
t2 = Sys.time()
out.joint <- fit.lasso(X, y, alpha=1, partition=0.6, penalty=penalty.fctr, 
                       n.iter=1000)
t <- Sys.time()
print((t - t2)) #2.987652 mins
save(out.joint, file = "Data/outJoint.Rda")






#####-------------------- OUTPUT: GLM --------------------#####
load("Data/outJoint.Rda")

evalGLM <- NULL
for(i in 1:length(out.joint)) {
  temp <- out.joint[[i]]$glm$evaluation
  temp$iter <- names(out.joint)[i]
  evalGLM <- rbind(evalGLM, temp)
}

ggplot(evalGLM[evalGLM$Type=="Test", ], aes(x=Balanced.Acc)) +
  geom_density() +
  geom_vline(aes(xintercept=mean(evalGLM$Balanced.Acc)), color="red", linetype="dashed", size=1) +
  annotate("text", x=mean(evalGLM$Balanced.Acc)-0.02, y=c(15), 
           label=paste0("Mean: ", round(mean(evalGLM$Balanced.Acc), 2)), 
           size=4.5, color = c("red"))+
  labs(x="Balanced Accuracy",
       y="Density") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))






#####-------------------- OUTPUT: LASSO --------------------#####
load("Data/outJoint.Rda")

evalLASSO <- NULL
for(i in 1:length(out.joint)) {
  temp <- out.joint[[i]]$lasso$evaluation
  temp$iter <- names(out.joint)[i]
  evalLASSO <- rbind(evalLASSO, temp)
}


# Microbiome frequency of selection
coef.lasso <- NULL
for(i in 1:length(out.joint)) {
  temp <- out.joint[[i]]$lasso$selected.var
  
  if(nrow(temp) == 0) {
    next
  } else {
    temp$iter <- names(out.joint)[i]
    coef.lasso <- rbind(coef.lasso, temp)
  }
  
  print(i)
}


# Remove those that is always 0
coef.lasso <- coef.lasso[abs(coef.lasso$coef) > 0, ]
coef.lasso <- coef.lasso[coef.lasso$var != "(Intercept)", ]

# Get the feature and genus name
load("Data/uBiome_taxtable.Rda")

ubiome <- unique(coef.lasso$var[grepl("svs", coef.lasso$var, ignore.case = TRUE)])
temp.ubiome <- coef.lasso[coef.lasso$var %in% ubiome, ]
temp.ubiome <- merge(temp.ubiome, taxTable, by.x = "var", by.y = "row.names", all.x = TRUE)
temp.ubiome$Genus <- as.factor(factor(temp.ubiome$Genus))

save(temp.ubiome, file = "Data/ubiomeSelect_joint.Rda")

order.x <- factor(unique(coef.lasso$Genus), levels = unique(coef.lasso$Genus))

give.n <- function(x, upper_limit = max(coef.lasso$coef, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}


d <- as.data.frame(table(temp.ubiome$Genus))
colnames(d) <- c("Genus", "freq")
d <- d[order(-d$freq), ]
d$Genus <- as.factor(factor(d$Genus, levels = unique(d$Genus[order(-d$freq)]), ordered=TRUE))

ggplot(data=d[1:20, ], aes(x=Genus, y=freq)) +
  geom_col() +
  # labs(x = "", y = "Selection proportion") +
  # scale_fill_manual(name = "> 0.5",
  #                   labels = c("No", "Yes"),
  #                   values = c("1" = "tomato3", "0" = "grey54")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle=90, hjust=1, size = 12))






#####-------------------- OUTPUT: PERMUTATION --------------------#####
load("Data/outJoint.Rda")

evalPERM <- NULL
for(i in 1:length(out.joint)) {
  temp <- out.joint[[i]]$permutation$evaluation
  temp$iter <- names(out.joint)[i]
  evalPERM <- rbind(evalPERM, temp)
}

ggplot(evalPERM[evalPERM$Type=="Test", ], aes(x=Balanced.Acc)) +
  geom_density() +
  geom_vline(aes(xintercept=mean(evalPERM$Balanced.Acc)), color="red", linetype="dashed", size=1) +
  annotate("text", x=mean(evalPERM$Balanced.Acc)-0.02, y=c(15), 
           label=paste0("Mean: ", round(mean(evalPERM$Balanced.Acc), 2)), 
           size=4.5, color = c("red"))+
  labs(x="Balanced Accuracy",
       y="Density") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))


# Microbiome frequency of selection
coef.perm <- NULL
for(i in 1:length(out.joint)) {
  temp <- out.joint[[i]]$permutation$selected.var
  
  if(nrow(temp) == 0) {
    next
  } else {
    temp$iter <- names(out.joint)[i]
    coef.perm <- rbind(coef.perm, temp)
  }
  
  print(i)
}


# Remove those that is always 0
coef.perm <- coef.perm[abs(coef.perm$coef) > 0, ]
coef.perm <- coef.perm[coef.perm$var != "(Intercept)", ]

# Get the feature and genus name
load("Data/uBiome_taxtable.Rda")

ubiome <- unique(coef.perm$var[grepl("svs", coef.perm$var, ignore.case = TRUE)])
temp.ubiome <- coef.perm[coef.perm$var %in% ubiome, ]
temp.ubiome <- merge(temp.ubiome, taxTable, by.x = "var", by.y = "row.names", all.x = TRUE)
temp.ubiome$Genus <- as.factor(factor(temp.ubiome$Genus))

save(temp.ubiome, file = "Data/ubiomeSelect_joint_permutation.Rda")

order.x <- factor(unique(coef.perm$Genus), levels = unique(coef.perm$Genus))

give.n <- function(x, upper_limit = max(coef.perm$coef, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}


d <- as.data.frame(table(temp.ubiome$Genus))
colnames(d) <- c("Genus", "freq")
d <- d[order(-d$freq), ]
d$Genus <- as.factor(factor(d$Genus, levels = unique(d$Genus[order(-d$freq)]), ordered=TRUE))

ggplot(data=d[1:20, ], aes(x=Genus, y=freq)) +
  geom_col() +
  # labs(x = "", y = "Selection proportion") +
  # scale_fill_manual(name = "> 0.5",
  #                   labels = c("No", "Yes"),
  #                   values = c("1" = "tomato3", "0" = "grey54")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle=90, hjust=1, size = 12))






#####-------------------- OUTPUT: LASSO VS GLMM --------------------#####
evalLASSO <- evalLASSO[evalLASSO$Type=="Test", ]
evalGLM <- evalGLM[evalGLM$Type=="Test", ]

dfLASSO <- evalLASSO[, names(evalLASSO) %in% c("Balanced.Acc", "iter")]
dfGLMM <- evalGLM[, names(evalGLM) %in% c("Balanced.Acc", "iter")]

names(dfLASSO)[1] <- "AccLASSO"
names(dfGLMM)[1] <- "AccGLMM"

lasso_glmm <- merge(dfLASSO, dfGLMM, by='iter')
lasso_glmm$gainLASSO <- (lasso_glmm$AccLASSO-lasso_glmm$AccGLMM)/lasso_glmm$AccLASSO


# Histogram of accuracy
library(data.table)
temp <- subset(lasso_glmm, select=-c(gainLASSO))
long_LassoGLMM <- melt(setDT(temp), id.vars=c("iter"), variable.name="Model")


library(dplyr)
library(magrittr)
stat_LassoGLMM <- long_LassoGLMM%>%
  group_by(Model)%>% 
  summarise(Median=median(value), Std=sd(value))

ggplot(long_LassoGLMM, aes(x=value, color=Model)) +
  geom_density(adjust=2) +
  scale_color_manual(values=c("blue", "red"),
                     breaks=c("AccLASSO", "AccGLMM")) +
  geom_vline(aes(xintercept=Median), stat_LassoGLMM, color=c("blue", "red"), linetype="dashed") +
  annotate("text", x=(stat_LassoGLMM$Median)+0.02, y=10,
           label=paste0("Median: ", round(stat_LassoGLMM$Median, 2)),
           size=4.5, color=c("blue", "red")) +
  labs(x="Accuracy",
       y="Density") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))


# Histogram of gaining of LASSO compared to GLMM
ggplot(lasso_glmm, aes(x=gainLASSO)) +
  geom_density(adjust=2) +
  geom_vline(aes(xintercept=mean(lasso_glmm$gainLASSO)), color="red", linetype="dashed", size=1) +
  annotate("text", x=median(lasso_glmm$gainLASSO)+0.02, y=c(15), 
           label=paste0("Median: ", round(median(lasso_glmm$gainLASSO), 2)), 
           size=4.5, color = c("red"))+
  labs(x="Accuracy Gained for LASSO",
       y="Density") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))






#####-------------------- OUTPUT: LASSO VS Permutation --------------------#####
evalLASSO <- evalLASSO[evalLASSO$Type=="Test", ]
evalPERM <- evalPERM[evalPERM$Type=="Test", ]

dfLASSO <- evalLASSO[, names(evalLASSO) %in% c("Balanced.Acc", "iter")]
dfPERM <- evalPERM[, names(evalPERM) %in% c("Balanced.Acc", "iter")]

names(dfLASSO)[1] <- "AccLASSO"
names(dfPERM)[1] <- "AccPERM"

lasso_perm <- merge(dfLASSO, dfPERM, by='iter')
lasso_perm$gainLASSO <- (lasso_perm$AccLASSO-lasso_perm$AccPERM)/lasso_perm$AccLASSO


# Histogram of accuracy
library(data.table)
temp <- subset(lasso_perm, select=-c(gainLASSO))
long_LassoPERM <- melt(setDT(temp), id.vars=c("iter"), variable.name="Model")


library(dplyr)
library(magrittr)
stat_LassoPERM <- long_LassoPERM%>%
  group_by(Model)%>% 
  summarise(Median=median(value), Std=sd(value))

ggplot(long_LassoPERM, aes(x=value, color=Model)) +
  geom_density(adjust=2) +
  scale_color_manual(values=c("blue", "red"),
                     breaks=c("AccLASSO", "AccPERM")) +
  geom_vline(aes(xintercept=Median), stat_LassoPERM, color=c("blue", "red"), linetype="dashed") +
  annotate("text", x=(stat_LassoPERM$Median)+0.02, y=10,
           label=paste0("Median: ", round(stat_LassoPERM$Median, 2)),
           size=4.5, color=c("blue", "red")) +
  labs(x="Accuracy",
       y="Density") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))


# Histogram of gaining of LASSO compared to Permutation
ggplot(lasso_perm, aes(x=gainLASSO)) +
  geom_density(adjust=2) +
  geom_vline(aes(xintercept=mean(lasso_perm$gainLASSO)), color="red", linetype="dashed", size=1) +
  annotate("text", x=median(lasso_perm$gainLASSO)+0.02, y=c(15), 
           label=paste0("Median: ", round(median(lasso_perm$gainLASSO), 2)), 
           size=4.5, color = c("red"))+
  labs(x="Accuracy Gained for LASSO",
       y="Density") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))






#####-------------------- OUTPUT: GLMM VS Permutation --------------------#####
evalGLM <- evalGLM[evalGLM$Type=="Test", ]
evalPERM <- evalPERM[evalPERM$Type=="Test", ]

dfGLMM <- evalGLM[, names(evalGLM) %in% c("Balanced.Acc", "iter")]
dfPERM <- evalPERM[, names(evalPERM) %in% c("Balanced.Acc", "iter")]

names(dfGLMM)[1] <- "AccGLMM"
names(dfPERM)[1] <- "AccPERM"

glmm_perm <- merge(dfGLMM, dfPERM, by='iter')
glmm_perm$gainGLMM <- (glmm_perm$AccGLMM-glmm_perm$AccPERM)/glmm_perm$AccGLMM


# Histogram of accuracy
library(data.table)
temp <- subset(glmm_perm, select=-c(gainGLMM))
long_glmmPERM <- melt(setDT(temp), id.vars=c("iter"), variable.name="Model")


library(dplyr)
library(magrittr)
stat_glmmPERM <- long_glmmPERM%>%
  group_by(Model)%>% 
  summarise(Median=median(value), Std=sd(value))

ggplot(long_glmmPERM, aes(x=value, color=Model)) +
  geom_density(adjust=2) +
  scale_color_manual(values=c("blue", "red"),
                     breaks=c("AccGLMM", "AccPERM")) +
  geom_vline(aes(xintercept=Median), stat_glmmPERM, color=c("blue", "red"), linetype="dashed") +
  annotate("text", x=(stat_glmmPERM$Median)+0.02, y=10,
           label=paste0("Median: ", round(stat_glmmPERM$Median, 2)),
           size=4.5, color=c("blue", "red")) +
  labs(x="Accuracy",
       y="Density") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))


# Histogram of gaining of LASSO compared to Permutation
ggplot(glmm_perm, aes(x=gainGLMM)) +
  geom_density(adjust=2) +
  geom_vline(aes(xintercept=mean(glmm_perm$gainGLMM)), color="red", linetype="dashed", size=1) +
  annotate("text", x=median(glmm_perm$gainGLMM)+0.02, y=c(15), 
           label=paste0("Median: ", round(median(glmm_perm$gainGLMM), 2)), 
           size=4.5, color = c("red"))+
  labs(x="Accuracy Gained for GLMM",
       y="Density") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))
