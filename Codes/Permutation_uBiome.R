####--------------------------------------------####
#### Permutation: 6 microbiome, 9695 microarray ####
####--------------------------------------------####

## Data preparation
load("~/02. PhD/CRC microbiome/drp_crcmicrobiome/Joint analysis/Data/rna_crcad.Rda")
load("~/02. PhD/CRC microbiome/drp_crcmicrobiome/Joint analysis/Data/otu_431_clr.Rda")
# load("~/02. PhD/CRC microbiome/drp_crcmicrobiome/Joint analysis/Data/uarraySelect_lasso_2.Rda")

genus <- get(load("~/02. PhD/CRC microbiome/drp_crcmicrobiome/Joint analysis/Data/selected genus.Rda"))
genus <- genus[order(-genus$freq), ]


# Get top 4 microbiome features
get.genus <- genus %>%
  filter(freq >= 500) %>%
  select(genus) %>%
  unlist() %>% unname()

otu.clr.4 <- otu.clr[, colnames(otu.clr) %in% as.character(get.genus)]


# Add microarray data
temp <- merge(otu.clr.4, rna.crcad, by="row.names")
rownames(temp) <- temp$Row.names
df.permut <- subset(temp, select=-c(Row.names))

save(df.permut, file="Data/df_permut_ubiome.Rda")





#####-------------------- START FROM HERE: CODES --------------------#####
setwd("/home/dputri/02. PhD/PhD_LungCancer")

load("~/02. PhD/CRC microbiome/drp_crcmicrobiome/Joint analysis/Data/group_crcad.Rda")
load("Data/df_permut_ubiome.Rda")
source("Codes/functions.R")


# Create the matrix for lasso input
y <- group
X <- data.matrix(df.permut)

uBiome <- colnames(X)[grepl("SVs", colnames(X), ignore.case = FALSE)]
uArray <- colnames(X)[grepl("_at", colnames(X), ignore.case = FALSE)]

penalty.fctr <- rep(0, ncol(X))
names(penalty.fctr) <- colnames(X)
penalty.fctr[uArray] <- 1


# LASSO
fit.lasso <- function(X, y, alpha=1, partition=0.6, penalty=NULL, 
                      n.iter=1000) {
  
  require(glmnet)
  require(caret)
  
  out <- vector("list", n.iter)
  
  
  for(i.k in 1:n.iter) {
    set.seed(i.k)
    options(warn = -1)
    
    
    # Create data partition
    ind <- createDataPartition(y, times = 1, p = partition, list = FALSE)
    
    
    # LOGISTIC REGRESSION
    dt.log <- data.frame(group=y[ind],
                         X[ind, uBiome])
    mod.glm <- glm(group ~ ., data=dt.log, family=binomial(link='logit'))
    out.glm <- summary.glm(mod.glm, x=as.data.frame(X), y=y, ind=ind) 
    
    
    # LASSO
    cv.lasso <- NULL
    while(is(cv.lasso, 'try-error') || is.null(cv.lasso)) {
      
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
      
    }
    
    
    mod.lasso <- glmnet(x=X[ind, ], y=y[ind], alpha=1, family="binomial",
                        penalty.factor=penalty, standardize=TRUE,
                        lambda = cv.lasso$lambda.min)
    
    out.lasso <- summary.lasso(mod.lasso, x=X, y=y, ind=ind)
    
    
    # PERMUTATED LASSO
    smpl <- sample(nrow(X))
    X.train.new <- data.frame(X[, uBiome], X[smpl, uArray])
    colnames(X.train.new)[5:ncol(X.train.new)] <- substring(colnames(X.train.new)[5:ncol(X.train.new)], 2)
    X.train.new <- data.matrix(X.train.new)
    
    cv.perm <- NULL
    while(is(cv.perm, 'try-error') || is.null(cv.perm)) {
      
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
      
    }
    
    mod.perm <- glmnet(x=X.train.new[ind, ], y=y[ind], alpha=1, family="binomial",
                       penalty.factor=penalty, standardize=TRUE,
                       lambda = cv.perm$lambda.min)
    
    out.perm <- summary.lasso(mod.perm, x=X.train.new, y=y, ind=ind)
    
    
    # Save the output
    out[[i.k]] <- list(out.glm, out.lasso, out.perm)
    names(out[[i.k]]) <- c("glm", "lasso", "permutation")
    names(out)[[i.k]] <- i.k
    
    print(i.k)
    
  }
  
  return(out)
  
}


# Standardized
set.seed(1)
t2 = Sys.time()
out.joint <- fit.lasso(X, y, alpha=1, partition=0.6, penalty=penalty.fctr, 
                       n.iter=500)
t <- Sys.time()
print((t - t2)) #1.126539 hours
save(out.joint, file = "Data/outPerm_uBiome.Rda")






#####-------------------- OUTPUT: GLM --------------------#####
load("Data/outPerm_uBiome.Rda")

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
load("Data/outPerm_uBiome.Rda")

evalLASSO <- NULL
for(i in 1:length(out.joint)) {
  temp <- out.joint[[i]]$lasso$evaluation
  temp$iter <- names(out.joint)[i]
  evalLASSO <- rbind(evalLASSO, temp)
}


# Microarray frequency of selection
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
load("~/02. PhD/CRC microbiome/drp_crcmicrobiome/Joint analysis/Data/esetRna_new.Rda") 

fdata <- as(featureData(esetRna), "data.frame")
temp <- subset(fdata, select = c(SYMBOL))

uArray <- colnames(X)[!grepl("SVs", colnames(X), ignore.case = TRUE)]

temp.uarray <- coef.lasso[coef.lasso$var %in% uArray, ]
temp.uarray$feat <- gsub("|`", "", temp.uarray$var)
temp.uarray <- merge(temp.uarray, temp, by.x = "var", by.y = "row.names", all.x = TRUE)
names(temp.uarray)[names(temp.uarray) == "SYMBOL"] <- "name"

save(temp.uarray, file = "Data/uarraySelect_lasso.Rda")

order.x <- factor(unique(coef.lasso$SYMBOL), levels = unique(coef.lasso$SYMBOL))

give.n <- function(x, upper_limit = max(coef.lasso$coef, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}


d <- as.data.frame(table(temp.uarray$name))
colnames(d) <- c("Genes", "Freq")
d <- d[order(-d$Freq), ]
d$Genes <- as.factor(factor(d$Genes, levels = unique(d$Genes[order(-d$Freq)]), ordered=TRUE))

ggplot(data=d[1:20, ], aes(x=Genes, y=Freq)) +
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


# Microarray frequency of selection
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

fdata <- as(featureData(esetRna), "data.frame")
temp <- subset(fdata, select = c(SYMBOL))

temp.uarray <- coef.perm[coef.perm$var %in% uArray, ]
temp.uarray <- merge(temp.uarray, temp, by.x = "var", by.y = "row.names", all.x = TRUE)
temp.uarray$SYMBOL <- as.factor(factor(temp.uarray$SYMBOL))

save(temp.uarray, file = "Data/uarraySelect_joint_permutation.Rda")

order.x <- factor(unique(coef.perm$SYMBOL), levels = unique(coef.perm$SYMBOL))

give.n <- function(x, upper_limit = max(coef.perm$coef, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}


d <- as.data.frame(table(temp.uarray$SYMBOL))
colnames(d) <- c("Genes", "freq")
d <- d[order(-d$freq), ]
d$Genes <- as.factor(factor(d$Genes, levels = unique(d$Genes[order(-d$freq)]), ordered=TRUE))

ggplot(data=d[1:20, ], aes(x=Genes, y=freq)) +
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
