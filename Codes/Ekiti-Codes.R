
load("Dataset2.RData")

#############################################################
#### Model Selection using stepwise logistic regression #####
#############################################################

library(MASS)

full.model<-glm(Indicator_Group~., family = "binomial", data=clinical)  #clinical contains response indicator and all 14 clinical factors
step<-stepAIC(full.model,trace=0, direction="backward") ##backward selection

final.model<-glm(Indicator_Group ~ Age + Packyears + BMI + COPD + Cardiac + Smoking + PDMT + HBP + Coagulation, family = "binomial",data = clinical)
Estimates<-summary(final.model)$coef ## I took Only significant ones


########################################################################
##                       UNIVARIATE aNALYSIS:                         ##
##  Fit logistic regression for clinical factors and one metabolite   ##
########################################################################

para<-pval<-numeric(102)
for (x in 1:ncol(Metabolic_data)){
model<-glm(metabolic$Indicator_Group~clinical$Age + clinical$Packyears + 
clinical$BMI + clinical$COPD + clinical$Cardiac + clinical$Smoking
+ clinical$Coagulation + Metabolic_data[[x]], family = "binomial")
para[x]<-summary(model)$coef[10,1]
pval[x]<-summary(model)$coef[10,4]
}

#### multiplicity correction 
library(multtest)

adjusted2 = mt.rawp2adjp(pval, "BH")
round(adjusted2$adjp, digits=6)
y<-log(adjusted2$adjp[,2])
plot(para,y,pch=19,col="blue", xlab="Coefficient Estimates",ylab="log(P-values)")
abline(a=NULL, b=NULL, h=log(0.05),col=2, v=NULL)

############################################################
##################    LASSO Regression    ##################

library(lattice)
library(glmnet)
library(coefplot)
library(caret)

trainY<-c(Response_data[,2])
trainx<-cbind(Clinical_data[,-c(2,6,7,8,11,12,13)],Metabolic_data)
trainX<-data.matrix(trainx, rownames.force = NA)

Lasso.fit<-cv.glmnet(x=trainX,y=trainY,family="binomial", type.measure="class",
alpha=1, nlambda=100)
plot(Lasso.fit)


Lasso.Coefs <- coef(Lasso.fit, s="lambda.min")
Est<-Lasso.Coefs[which(Lasso.Coefs != 0 )] 
Lasso.Coefs@Dimnames[[1]][which(Lasso.Coefs != 0 )]
round(Est, digits=4)

#######################################################
#####     LASSO Regression with no penalty       ######
#####           for clinical factors             ######
          

penalty <- c(rep(0, 8), rep(1, 102)) #penalty factor

Lasso.fit2<-cv.glmnet(x=trainX,y=as.factor(trainY),family="binomial", type.measure="class",
alpha=1, nlambda=100,penalty.factor=penalty)

Lasso.Coefs2 <- coef(Lasso.fit2, s="lambda.1se")
Est2<-Lasso.Coefs[which(Lasso.Coefs2 != 0 ) ] 
Lasso.Coefs2@Dimnames[[1]][which(Lasso.Coefs2 != 0 ) ]
round(Est2, digits=4)
plot(Lasso.fit2)
###### Classification 
p2<-predict(Lasso.fit2,trainX,s="lambda.min", type="class")
confusionMatrix(as.factor(p2),as.factor(trainY))

##################################################################
#### Testing the additional predictive value of metabolic data####

           
         ####        Global Test I      ####

library(mboost)
library(stabs)
library(survival)
library(globalboosttest)

Mdata<-data.matrix(Metabolic_data, rownames.force = NA)
Cdata<-data.matrix(Clinical_data, rownames.force = NA)

Cdata<-Clinical_data[,-c(2,6,7,8,11,12,13)]
Test<-globalboosttest(Mdata,trainY,Z=Cdata,nperm=1000,mstop=1000,mstopAIC=FALSE,pvalueonly=F,plot=FALSE)

Test$pvalue


         ####        Global Test II     ####


library(globaltest)


GT<-gt(Indicator_Group ~ Age + Packyears + BMI + COPD + Cardiac + Smoking + Coagulation, ~., data = all.final)
summary(GT)



        ####        Likelihood Ratio Test     ####

full.model<-glm(Indicator_Group ~., family = "binomial",data = all.final)

reduced.model<-glm(Indicator_Group ~ Age + Packyears + BMI + COPD + Cardiac +
Smoking + Coagulation, family = "binomial",data = all.final)

LRT<-2*(logLik(full.model)-logLik(reduced.model)) #411.3724 
p.value<-1-pchisq(LRT,102)


##############################################################################
########                                                          ############
########                   Final Model Combined                   ############
##############################################################################

ALL.train<-all.final
library(glmnet)

ALL.train[,7]<-factor(ALL.train[,7],c("Never","Active","Stopped"))
table(ALL.train[,7])

y=ALL.train[,1] ##### THE RESPONSE  TRAIN ####
levels(y)<-c("0","1")
ALL.train2<-data.frame(ALL.train) #ALL.train2=ALL.train
h2=model.matrix(~.,ALL.train2[,-1])
p.fac2 = rep(1,dim(h2[,-1])[2])
k=ncol(h2[,-1])-ncol(ALL.train[,c(8:110)])
p.fac2[c(1:k)] = 0
X<-h2[,-1]  

##########################  FUNCTION ################################################
Fit.Lasso.kCV<-function(X,y,ind.train=1:nrow(X),ind.test=1:nrow(X),alpha.value=1,penality=NULL,indipend.data){
	options(warn=-1)
	if(is.null(penality)){
		fit21 <- try(cv.glmnet(X[ind.train,],y[ind.train],alpha = alpha.value,standardize = F,family="binomial",type.measure="class",nfold=3),silent=T)
		pred.fit21 <- try(predict(fit21,X[ind.test,],type="class"),silent=T)
		pred.fit213 <-try(predict(fit21,X[ind.train,],type="class"),silent=T)
            Coefficients <-try(coef(fit21, s="lambda.min"),silent=T)
		
		MCE.train<-try(1-sum(diag(table(y[ind.train],pred.fit213)))/length(y[ind.train]),silent=T)
		MCE<-try(1-sum(diag(table(y[ind.test],pred.fit21)))/length(y[ind.test]),silent=T)
		
		return(list(MCEi=MCE,MCE.train=MCE.train,pred.train=pred.fit213,pred.test=pred.fit21,train=ind.train,test=ind.test,Coefficients=Coefficients))
	}
	else {
		
		fit21 <- try(cv.glmnet(X[ind.train,],y[ind.train],intercept=T,alpha = alpha.value,standardize = F,family="binomial",type.measure="class", penalty.factor =penality,nfold=3),silent=T)
		pred.fit21 <- try(predict(fit21,X[ind.test,],type="class"),silent=T)
		pred.fit213 <- try(predict(fit21,X[ind.train,],type="class"),silent=T)
            Coefficients <-try(coef(fit21, s="lambda.min"),silent=T)
           
		MCE.train<-try(1-sum(diag(table(y[ind.train],pred.fit213)))/length(y[ind.train]),silent=T)
		MCE<-try(1-sum(diag(table(y[ind.test],pred.fit21)))/length(y[ind.test]),silent=T)
		
		return(list(MCEi=MCE,MCE.train=MCE.train,pred.train=pred.fit213,pred.test=pred.fit21,train=ind.train,test=ind.test,Coefficients=Coefficients))	
	}
}


####################### The LASSO Analysis ##############################
n.cv<-1000
MCE.test<-rep(NA,n.cv)
MCE.train<-rep(NA,n.cv)
Pred.train<-matrix(NA,357,n.cv)
Pred.test<-matrix(NA,179,n.cv)
Coeff<-matrix(NA,111,n.cv)
train<-matrix(NA,357,n.cv)
test<-matrix(NA,179,n.cv)

set.seed(110) 
for (j in 1:n.cv){ 
	ind.train <-as.vector(sort(sample(1:nrow(X),floor(2*nrow(X)/3),replace=F) ) )
	ind.test <- as.vector(c(1:nrow(X))[-ind.train])
	res.cv<-try(Fit.Lasso.kCV(X=X,y=y,ind.train=ind.train ,ind.test=ind.test,alpha.value=1, penality=p.fac2),silent=T)
	MCE.test[j]<-try(res.cv$MCEi,silent=T)
	MCE.train[j]<-try(res.cv$MCE.train,silent=T)
	Pred.train[,j]<-try(as.numeric(res.cv$pred.train),silent=T)
	Pred.test[,j]<-try(as.numeric(res.cv$pred.test),silent=T)
	train[,j]<-ind.train
	test[,j]<-ind.test
	Coeff[,j]<-(try((res.cv$Coefficients)[1:111], silent=T))
	print(j)
}
save(MCE.train,MCE.test,train,test,Pred.test,Pred.train,Coeff,file="3foldCV_ALL_Data.Rdata")

npv.test<-ppv.test<-sen.test<-spe.test<-rep(NA,n.cv)

for (j in 1:n.cv){ 
	tab.test<-table(y[test[,j]],Pred.test[,j])
	
	npv.test[j]<-tab.test[1,1]/sum(tab.test[1,1],tab.test[2,1])
	ppv.test[j]<-tab.test[2,2]/sum(tab.test[1,2],tab.test[2,2])
	sen.test[j]<-tab.test[2,2]/sum(tab.test[2,1],tab.test[2,2])
	spe.test[j]<-tab.test[1,1]/sum(tab.test[1,2],tab.test[1,1])
	print(j)	
}
save(MCE.train,MCE.test,train,test,npv.test,ppv.test,sen.test,spe.test,Coeff,file="3foldCV_ALL_Data_FINAL.Rdata")


##############################################################################
########                                                          ############
########                Final Model Clinical Data                 ############
##############################################################################
ALL.train<-all.final[,1:8]

ALL.train[,7]<-factor(ALL.train[,7],c("Never","Active","Stopped"))
table(ALL.train[,7])

y=ALL.train[,1] ##### THE RESPONSE  TRAIN ####
levels(y)<-c("0","1")
ALL.train2<-data.frame(ALL.train) #ALL.train2=ALL.train
h2=model.matrix(~.,ALL.train2[,-1])
p.fac2=c(0,0,0,0,0,0,0,1) ## I penalized "Coagulation" in both models instead of removing it
X<-h2[,-1]  

##########################  FUNCTION ################################################
Fit.Lasso.kCV<-function(X,y,ind.train=1:nrow(X),ind.test=1:nrow(X),alpha.value=1,penality=NULL,indipend.data){
	options(warn=-1)
	if(is.null(penality)){
		fit21 <- try(cv.glmnet(X[ind.train,],y[ind.train],alpha = alpha.value,standardize = F,family="binomial",type.measure="class",nfold=3),silent=T)
		pred.fit21 <- try(predict(fit21,X[ind.test,],type="class"),silent=T)
		pred.fit213 <-try(predict(fit21,X[ind.train,],type="class"),silent=T)
            Coefficients <-try(coef(fit21, s="lambda.min"),silent=T)
		
		MCE.train<-try(1-sum(diag(table(y[ind.train],pred.fit213)))/length(y[ind.train]),silent=T)
		MCE<-try(1-sum(diag(table(y[ind.test],pred.fit21)))/length(y[ind.test]),silent=T)
		
		return(list(MCEi=MCE,MCE.train=MCE.train,pred.train=pred.fit213,pred.test=pred.fit21,train=ind.train,test=ind.test,Coefficients=Coefficients))
	}
	else {
		
		fit21 <- try(cv.glmnet(X[ind.train,],y[ind.train],intercept=T,alpha = alpha.value,standardize = F,family="binomial",type.measure="class", penalty.factor =penality,nfold=3),silent=T)
		pred.fit21 <- try(predict(fit21,X[ind.test,],type="class"),silent=T)
		pred.fit213 <- try(predict(fit21,X[ind.train,],type="class"),silent=T)
            Coefficients <-try(coef(fit21, s="lambda.min"),silent=T)
           
		MCE.train<-try(1-sum(diag(table(y[ind.train],pred.fit213)))/length(y[ind.train]),silent=T)
		MCE<-try(1-sum(diag(table(y[ind.test],pred.fit21)))/length(y[ind.test]),silent=T)
		
		return(list(MCEi=MCE,MCE.train=MCE.train,pred.train=pred.fit213,pred.test=pred.fit21,train=ind.train,test=ind.test,Coefficients=Coefficients))	
	}
}


############################################## The LASSO Analysis ###################################################
n.cv<-1000
MCE.test<-rep(NA,n.cv)
MCE.train<-rep(NA,n.cv)
Pred.train<-matrix(NA,357,n.cv)
Pred.test<-matrix(NA,179,n.cv)
train<-matrix(NA,357,n.cv)
test<-matrix(NA,179,n.cv)

set.seed(110) 
for (j in 1:n.cv){ 
	ind.train <-as.vector(sort(sample(1:nrow(X),floor(2*nrow(X)/3),replace=F) ) )
	ind.test <- as.vector(c(1:nrow(X))[-ind.train])
	res.cv<-try(Fit.Lasso.kCV(X=X,y=y,ind.train=ind.train ,ind.test=ind.test,alpha.value=1, penality=p.fac2),silent=T)
	MCE.test[j]<-try(res.cv$MCEi,silent=T)
	MCE.train[j]<-try(res.cv$MCE.train,silent=T)
	Pred.train[,j]<-try(as.numeric(res.cv$pred.train),silent=T)
	Pred.test[,j]<-try(as.numeric(res.cv$pred.test),silent=T)
	train[,j]<-ind.train
	test[,j]<-ind.test
	#Coeff[,j]<-(try((res.cv$Coefficients)[1:111], silent=T))
	print(j)
}
save(MCE.train,MCE.test,train,test,Pred.test,Pred.train,file="3foldCV_CLINICAL_Data.Rdata")


npv.test<-ppv.test<-sen.test<-spe.test<-rep(NA,n.cv)

for (j in 1:n.cv){ 
	tab.test<-table(y[test[,j]],Pred.test[,j])
	
	npv.test[j]<-tab.test[1,1]/sum(tab.test[1,1],tab.test[2,1])
	ppv.test[j]<-tab.test[2,2]/sum(tab.test[1,2],tab.test[2,2])
	sen.test[j]<-tab.test[2,2]/sum(tab.test[2,1],tab.test[2,2])
	spe.test[j]<-tab.test[1,1]/sum(tab.test[1,2],tab.test[1,1])
	print(j)	
}
save(MCE.train,MCE.test,train,test,npv.test,ppv.test,sen.test,spe.test,file="3foldCV_CLINICAL_Data_FINAL.Rdata")





