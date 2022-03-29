####--------------------------------------------####
#### Permutation: 6 microarray, 431 microbiome ####
####--------------------------------------------####
load("Data/group_crcad.Rda")
load("Data/df_lasso_uArray.Rda")
source("Codes/functions.R")



## Load libraries
library(parallel)
library(ggplot2)
library(glmnet)





#####-------------------- CODES --------------------#####
dt <- as.matrix.data.frame(df.uArray[, names(df.uArray)!="group"])

uBiome <- colnames(dt)[grepl("SVs", colnames(dt), ignore.case = FALSE)]
uArray <- colnames(dt)[grepl("_at", colnames(dt), ignore.case = FALSE)]

penalty.fctr <- rep(0, ncol(dt))
names(penalty.fctr) <- colnames(dt)
penalty.fctr[uBiome] <- 1


# Create a cluster of cores
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl=cl, 
              varlist=c("dt", "group", "penalty.fctr", 
                        "uArray", "uBiome",
                        "summary.lasso", "summary.glm"))


# Main algorithm
system.time({
  out <- parLapply(cl=cl,
                   X=1:1000,
                   fun=function(i) {
                     require(glmnet)
                     require(caret)
                     
                     # sample 2/3 of the subjects at random
                     id_keep <- sample(x=rownames(dt),
                                       size=floor(2*nrow(dt)/3),
                                       replace=FALSE)
                     
                     # LOGISTIC REGRESSION
                     dt.log <- data.frame(group=group[id_keep],
                                          dt[id_keep, uArray])
                     names(dt.log)[2:ncol(dt.log)] <- substring(names(dt.log)[2:ncol(dt.log)], 2)
                     mod.glm <- glm(group ~ ., data=dt.log, family=binomial(link='logit'))
                     out.glm <- summary.glm(mod.glm, x=as.data.frame(dt), y=group, ind=id_keep) 
                     
                     
                     # LASSO
                     cv.lasso <- try(cv.glmnet(x=dt[id_keep, ], y=group[id_keep], alpha=1,
                                               standardize=TRUE, family="binomial",
                                               penalty.factor=penalty.fctr, 
                                               nfolds=5, type.measure="class"),
                                     silent=TRUE)
                     
                     m1 <- glmnet(x=dt[id_keep, ], y=group[id_keep], alpha=1, family="binomial",
                                  penalty.factor=penalty.fctr, standardize=TRUE,
                                  lambda=cv.lasso$lambda.min)
                     
                     out.lasso <- summary.lasso(m1, x=dt, y=group, ind=id_keep)
                     
                     
                     # PERMUTATED LASSO
                     smpl <- sample(nrow(dt))
                     X.train.new <- cbind(dt[, uArray], dt[smpl, uBiome])
                     # colnames(X.train.new)[1:6] <- substring(colnames(X.train.new)[1:6], 2)
                     X.train.new <- data.matrix(X.train.new)
                     
                     cv.perm <- try(cv.glmnet(x=X.train.new[id_keep, ], y=group[id_keep], alpha=1,
                                              standardize=TRUE, family="binomial",
                                              penalty.factor=penalty.fctr, 
                                              nfolds=5, type.measure="class"),
                                    silent = TRUE)
                     
                     mod.perm <- glmnet(x=X.train.new[id_keep, ], y=group[id_keep], alpha=1, family="binomial",
                                        penalty.factor=penalty.fctr, standardize=TRUE,
                                        lambda=cv.perm$lambda.min)
                     
                     out.perm <- summary.lasso(mod.perm, x=X.train.new, y=group, ind=id_keep)
                     
                     
                     # Save the output
                     out <- list(out.glm, out.lasso, out.perm)
                     names(out) <- c("glm", "lasso", "permutation")
                     
                     return(out)
                   }
  )
})


# Stop the cluster
stopCluster(cl)
gc() # Running time: 26.090s

save(out, file="Data/perm_uArray.Rda")







#####-------------------- OUTPUT: GLM --------------------#####
load("Data/perm_uArray.Rda")

library(tibble)
library(tidyr)
library(magrittr)
library(dplyr)


# Evaluation
tmp <- lapply(out, function(a) { #out[[1]]
  x <- lapply(a, function(b) { #out[[1]]$glm
    b$evaluation
  })
  
  y <- do.call("rbind", x) %>%
    rownames_to_column(var="model")
})

eval <- do.call("rbind", tmp) %>%
  mutate(model=do.call("rbind", strsplit(model, "[.]"))[, 1])

save(eval, file="Data/perf eval_permut uArray.Rda")


# Plot
ggplot(eval[eval$Type=="Test" & eval$model=="glm", ], aes(x=Balanced.Acc)) +
  geom_density() +
  geom_vline(aes(xintercept=mean(eval$Balanced.Acc)), color="red", linetype="dashed", size=1) +
  annotate("text", x=mean(eval$Balanced.Acc)-0.02, y=c(15), 
           label=paste0("Mean: ", round(mean(eval$Balanced.Acc), 2)), 
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
load("Data/perf eval_permut uArray.Rda")
load("Data/perm_uArray.Rda")


## Evaluation
ggplot(eval[eval$Type=="Test" & eval$model=="lasso", ], aes(x=Balanced.Acc)) +
  geom_density() +
  geom_vline(aes(xintercept=mean(eval$Balanced.Acc)), color="red", linetype="dashed", size=1) +
  annotate("text", x=mean(eval$Balanced.Acc)-0.02, y=c(15), 
           label=paste0("Mean: ", round(mean(eval$Balanced.Acc), 2)), 
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



## Feature selection
tmp <- lapply(out, function(a) { #out[[1]]
  a$lasso$selected.var
})

coef <- do.call("rbind", tmp) %>%
  filter(abs(coef) > 0,
         var!="(Intercept)")
rownames(coef) <- NULL

coef <- coef %>%
  select(-coef) %>%
  table %>%
  data.frame() %>%
  set_colnames(c("Feat", "Freq")) 



## Get the genes and genus name
library(Biobase)


# uArray
load("Data/esetRna_new.Rda") 

fdata <- as(featureData(esetRna), "data.frame")
temp <- subset(fdata, select = c(SYMBOL))

coef <- merge(coef, temp, by.x="Feat", by.y="row.names", all.x=TRUE) %>%
  rename(name="SYMBOL")


# uBiome
load("Data/uBiome_taxtable.Rda")

coef <- merge(coef, taxTable, by.x="Feat", by.y="row.names", all.x=TRUE) %>%
  mutate(name=ifelse(is.na(name), as.character(Genus), name)) %>%
  select(-Genus) %>%
  mutate(prop=Freq/1000)

save(coef, file = "Data/coef_permut uArray.Rda")


d <- coef %>%
  filter(Feat %in% Feat[grepl("SVs", Feat)])
d <- d[order(-d$Freq), ]


ggplot(data=d[1:20, ], aes(x=reorder(name, -Freq), y=Freq)) +
  geom_col() +
  labs(x = "", y = "Frequency of selection") +
  # scale_fill_manual(name = "> 0.5",
  #                   labels = c("No", "Yes"),
  #                   values = c("1" = "tomato3", "0" = "grey54")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle=90, hjust=1, size = 12))






#####-------------------- OUTPUT: PERMUTATION --------------------#####
load("Data/perm_uArray.Rda")
load("Data/perf eval_permut uArray.Rda")


## Evaluation
ggplot(eval[eval$Type=="Test" & eval$model=="permutation", ], aes(x=Balanced.Acc)) +
  geom_density() +
  geom_vline(aes(xintercept=mean(eval$Balanced.Acc)), color="red", linetype="dashed", size=1) +
  annotate("text", x=mean(eval$Balanced.Acc)-0.02, y=c(15), 
           label=paste0("Mean: ", round(mean(eval$Balanced.Acc), 2)), 
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



## Feature selection
tmp <- lapply(out, function(a) { #out[[1]]
  a$permutation$selected.var
})

coef <- do.call("rbind", tmp) %>%
  filter(abs(coef) > 0,
         var!="(Intercept)")
rownames(coef) <- NULL

coef <- coef %>%
  select(-coef) %>%
  table %>%
  data.frame() %>%
  set_colnames(c("Feat", "Freq")) 



## Get the genes and genus name
library(Biobase)


# uArray
load("Data/esetRna_new.Rda") 

fdata <- as(featureData(esetRna), "data.frame")
temp <- subset(fdata, select = c(SYMBOL))

coef <- merge(coef, temp, by.x="Feat", by.y="row.names", all.x=TRUE) %>%
  rename(name="SYMBOL")


# uBiome
load("Data/uBiome_taxtable.Rda")

coef <- merge(coef, taxTable, by.x="Feat", by.y="row.names", all.x=TRUE) %>%
  mutate(name=ifelse(is.na(name), as.character(Genus), name)) %>%
  select(-Genus) %>%
  mutate(prop=Freq/1000)

save(coef, file = "Data/coef permut_permut uArray.Rda")


d <- coef %>%
  filter(Feat %in% Feat[grepl("SVs", Feat)])
d <- d[order(-d$Freq), ]


ggplot(data=d[1:20, ], aes(x=reorder(name, -Freq), y=Freq)) +
  geom_col() +
  labs(x = "", y = "Frequency of selection") +
  # scale_fill_manual(name = "> 0.5",
  #                   labels = c("No", "Yes"),
  #                   values = c("1" = "tomato3", "0" = "grey54")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle=90, hjust=1, size = 12))






#####-------------------- OUTPUT: LASSO VS GLMM --------------------#####
evalLASSO <- eval[eval$Type=="Test" & eval$model=="lasso", ] %>%
  select(Balanced.Acc) %>%
  rename(AccLASSO="Balanced.Acc")
evalGLM <- eval[eval$Type=="Test" & eval$model=="glm", ] %>%
  select(Balanced.Acc) %>%
  rename(AccGLMM="Balanced.Acc")

lasso_glmm <- cbind(evalLASSO, evalGLM)
lasso_glmm$gainLASSO <- (lasso_glmm$AccLASSO-lasso_glmm$AccGLMM)/lasso_glmm$AccLASSO


# Histogram of accuracy
library(data.table)
library(reshape2)

temp <- subset(lasso_glmm, select=-c(gainLASSO))
temp$iter <- seq(1, nrow(temp))
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
  geom_vline(aes(xintercept=median(lasso_glmm$gainLASSO)), color="red", linetype="dashed", size=1) +
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
evalLASSO <- eval[eval$Type=="Test" & eval$model=="lasso", ] %>%
  select(Balanced.Acc) %>%
  rename(AccLASSO="Balanced.Acc")
evalPERM <- eval[eval$Type=="Test" & eval$model=="permutation", ] %>%
  select(Balanced.Acc) %>%
  rename(AccPERM="Balanced.Acc")

lasso_perm <- cbind(evalLASSO, evalPERM)
lasso_perm$gainLASSO <- (lasso_perm$AccLASSO-lasso_perm$AccPERM)/lasso_perm$AccLASSO


# Histogram of accuracy
library(data.table)
library(reshape2)

temp <- subset(lasso_perm, select=-c(gainLASSO))
temp$iter <- seq(1, nrow(temp))
long_LassoPERM <- melt(setDT(temp), id.vars=c("iter"), variable.name="Model")


library(dplyr)
library(magrittr)
stat_LassoPERM <- long_LassoPERM%>%
  group_by(Model)%>% 
  summarise(Median=median(value), Std=sd(value))

ggplot(long_LassoGLMM, aes(x=value, color=Model)) +
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


# Histogram of gaining of LASSO compared to GLMM
ggplot(lasso_perm, aes(x=gainLASSO)) +
  geom_density(adjust=2) +
  geom_vline(aes(xintercept=median(lasso_perm$gainLASSO)), color="red", linetype="dashed", size=1) +
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
evalGLM <- eval[eval$Type=="Test" & eval$model=="glm", ] %>%
  select(Balanced.Acc) %>%
  rename(AccGLM="Balanced.Acc")
evalPERM <- eval[eval$Type=="Test" & eval$model=="permutation", ] %>%
  select(Balanced.Acc) %>%
  rename(AccPERM="Balanced.Acc")

GLM_perm <- cbind(evalGLM, evalPERM)
GLM_perm$gainGLM <- (GLM_perm$AccGLM-GLM_perm$AccPERM)/GLM_perm$AccGLM


# Histogram of accuracy
library(data.table)
library(reshape2)

temp <- subset(GLM_perm, select=-c(gainGLM))
temp$iter <- seq(1, nrow(temp))
long_GLMPERM <- melt(setDT(temp), id.vars=c("iter"), variable.name="Model")


library(dplyr)
library(magrittr)
stat_GLMPERM <- long_GLMPERM%>%
  group_by(Model)%>% 
  summarise(Median=median(value), Std=sd(value))

ggplot(long_GLMPERM, aes(x=value, color=Model)) +
  geom_density(adjust=2) +
  scale_color_manual(values=c("blue", "red"),
                     breaks=c("AccGLM", "AccPERM")) +
  geom_vline(aes(xintercept=Median), stat_GLMPERM, color=c("blue", "red"), linetype="dashed") +
  annotate("text", x=(stat_GLMPERM$Median)+0.02, y=10,
           label=paste0("Median: ", round(stat_GLMPERM$Median, 2)),
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


# Histogram of gaining of GLM compared to PERMUTATION
ggplot(GLM_perm, aes(x=gainGLM)) +
  geom_density(adjust=2) +
  geom_vline(aes(xintercept=median(GLM_perm$gainGLM)), color="red", linetype="dashed", size=1) +
  annotate("text", x=median(GLM_perm$gainGLM)+0.02, y=c(15), 
           label=paste0("Median: ", round(median(GLM_perm$gainGLM), 2)), 
           size=4.5, color = c("red"))+
  labs(x="Accuracy Gained for GLM",
       y="Density") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))
