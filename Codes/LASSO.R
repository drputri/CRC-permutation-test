###############  ALL DATA ##########
library(glmnet)

load(file="ALLDATA.Rdata")

clin<-c("Group","Age","Packyears","BMI","COPD","Cardiac","Smoking","Coagulation") #CHOSEN USING BACKWARD SELECTION AND KEEPING SIGNIFICANT ONES

ALL.train[,11]<-factor(ALL.train[,11],c("Never","Active","Stopped"))
table(ALL.train[,11])
ALL.test[,11]<-factor(ALL.test[,11],c("Never","Active","Stopped"))
table(ALL.test[,11])

y=ALL.train[,1] ##### THE RESPONSE  TRAIN ####
levels(y)<-c("0","1")
ALL.train2<-data.frame(ALL.train[,clin],ALL.train[,c(16:117)])
h2=model.matrix(~.,ALL.train2[,-1])
p.fac2 = rep(1,dim(h2[,-1])[2])
k=ncol(h2[,-1])-ncol(ALL.train[,c(16:117)])
p.fac2[c(1:k)] = 0
X<-h2[,-1]  
y2=ALL.test[,1] ##### THE RESPONSE  TRAIN ####
levels(y2)<-c("0","1")
ALL.test2<-data.frame(ALL.test[,clin],ALL.test[,c(16:117)])

h3=model.matrix(~.,ALL.test2)
indip<-h3[,-1]
#----------------------  3-FOLD CV for LASSO ---------------------

############################### MODIFIED FUNCTION #####################
### INDEPENDENT.DATA CONTAINS THE TEST DATA [FIRST COLUMN IS Y]
Fit.Lasso.kCV<-function(X,y,ind.train=1:nrow(X),ind.test=1:nrow(X),alpha.value=1,penality=NULL,indipend.data){
	options(warn=-1)
	if(is.null(penality)){
		fit21 <- try(cv.glmnet(X[ind.train,],y[ind.train],alpha = alpha.value,standardize = F,family="binomial",type.measure="class",nfold=3),silent=T)
		pred.fit21 <- try(predict(fit21,X[ind.test,],type="class"),silent=T)
		pred.fit212<- try(predict(fit21,indipend.data[,-1],type="class"),silent=T)
		pred.fit213 <-try(predict(fit21,X[ind.train,],type="class"),silent=T)
		
		MCEtrain<-try(1-sum(diag(table(y[ind.train],pred.fit213)))/length(y[ind.train]),silent=T)
		MCE<-try(1-sum(diag(table(y[ind.test],pred.fit21)))/length(y[ind.test]),silent=T)
		MCE.test<-try(1-sum(diag(table(indipend.data[,1],pred.fit212)))/length(indipend.data[,1]),silent=T)
		
		
		
		return(list(MSEi=MCE,MCE.test=MCE.test,MCEtrain=MCEtrain,pred.test=pred.fit212,pred.train=pred.fit213,pred.data=pred.fit21,train=ind.train,test=ind.test))
	}
	else {
		
		fit21 <- try(cv.glmnet(X[ind.train,],y[ind.train],intercept=T,alpha = alpha.value,standardize = F,family="binomial",type.measure="class",penalty.factor =penality,nfold=3),silent=T)
		pred.fit21 <- try(predict(fit21,X[ind.test,],type="class"),silent=T)
		pred.fit212 <- try(predict(fit21,indipend.data[,-1],type="class"),silent=T)
		pred.fit213 <- try(predict(fit21,X[ind.train,],type="class"),silent=T)
		
		MCEtrain<-try(1-sum(diag(table(y[ind.train],pred.fit213)))/length(y[ind.train]),silent=T)
		MCE<-try(1-sum(diag(table(y[ind.test],pred.fit21)))/length(y[ind.test]),silent=T)
		MCE.test<-try(1-sum(diag(table(indipend.data[,1],pred.fit212)))/length(indipend.data[,1]),silent=T)
		
		
		
		return(list(MSEi=MCE,MCE.test=MCE.test,MCEtrain=MCEtrain,pred.test=pred.fit212,pred.train=pred.fit213,pred.data=pred.fit21,train=ind.train,test=ind.test))	
	}
}





######################## LASOOOOO########################

n.cv<-1000 # no.of cross validations
MSE.data<-rep(NA,n.cv)
MSE.test<-rep(NA,n.cv)
MCEtrain<-rep(NA,n.cv)
Pred.test<-matrix(NA,168,n.cv)
Pred.train<-matrix(NA,357,n.cv)
Pred.data<-matrix(NA,179,n.cv)

train<-matrix(NA,357,n.cv)
test<-matrix(NA,179,n.cv)

set.seed(110) 

for (j in 1:n.cv){ 
	ind.train <-as.vector(sort(sample(1:nrow(X),floor(2*nrow(X)/3),replace=F) ) )
	ind.test <- as.vector(c(1:nrow(X))[-ind.train])
	res.cv<-try(Fit.Lasso.kCV(X=X,y=y,ind.train=ind.train ,ind.test=ind.test,alpha.value=1,penality=p.fac2,indipend.data=indip),silent=T)
	#res.cv<-Fit.Lasso.kCV(X=X,y=y,ind.train=ind.train ,ind.test=ind.test,alpha.value=1,indipend.data=indip)
	MSE.data[j]<-try(res.cv$MSEi,silent=T)
	MSE.test[j]<-try(res.cv$MCE.test,silent=T)
	MCEtrain[j]<-try(res.cv$MCEtrain,silent=T)
	Pred.test[,j]<-try(as.numeric(res.cv$pred.test),silent=T)
	Pred.train[,j]<-try(as.numeric(res.cv$pred.train),silent=T)
	Pred.data[,j]<-try(as.numeric(res.cv$pred.data),silent=T)
	train[,j]<-ind.train
	test[,j]<-ind.test
	
	print(j)
	
}




save(MCEtrain,MSE.test,MSE.data,train,test,Pred.data,Pred.train,Pred.test,file="3foldCV_LASSO(ALLDATA noPneality).Rdata")


npv.test<-ppv.test<-sen.test<-spe.test<-rep(NA,n.cv)
npv.data<-ppv.data<-sen.data<-spe.data<-rep(NA,n.cv)


for (j in 1:n.cv){ 
	
	
	tab.data<-table(y[test[,j]],Pred.data[,j])
	
	npv.data[j]<-tab.data[1,1]/sum(tab.data[1,1],tab.data[2,1])
	ppv.data[j]<-tab.data[2,2]/sum(tab.data[1,2],tab.data[2,2])
	sen.data[j]<-tab.data[2,2]/sum(tab.data[2,1],tab.data[2,2])
	spe.data[j]<-tab.data[1,1]/sum(tab.data[1,2],tab.data[1,1])
	
	tab.test<-table(y2,Pred.test[,j])
	
	npv.test[j]<-tab.test[1,1]/sum(tab.test[1,1],tab.test[2,1])
	ppv.test[j]<-tab.test[2,2]/sum(tab.test[1,2],tab.test[2,2])
	sen.test[j]<-tab.test[2,2]/sum(tab.test[2,1],tab.test[2,2])
	spe.test[j]<-tab.test[1,1]/sum(tab.test[1,2],tab.test[1,1])
	
	print(j)
	
}



save(MCEtrain,MSE.test,MSE.data,train,test,npv.test,ppv.test,sen.test,spe.test,npv.data,
     ppv.data,sen.data,spe.data,file="3foldCV_LASSO(ALLDATA noPneality)FINAL.Rdata")









######################## LASSO   PERMUTATION  ########################

n.cv<-1000#00 # no.of cross validations
MSE.data<-rep(NA,n.cv)
MSE.test<-rep(NA,n.cv)
MCEtrain<-rep(NA,n.cv)
Pred.test<-matrix(NA,168,n.cv)
Pred.train<-matrix(NA,357,n.cv)
Pred.data<-matrix(NA,179,n.cv)

train<-matrix(NA,357,n.cv)
test<-matrix(NA,179,n.cv)

set.seed(110) 

for (j in 1:n.cv){ 
	ind.train <-as.vector(sort(sample(1:nrow(X),floor(2*nrow(X)/3),replace=F) ) )
	ind.test <- as.vector(c(1:nrow(X))[-ind.train])
	n<-nrow(X)
	set.seed(j)
	sampi <- sample(n)
	ALL.train2<-data.frame(ALL.train[,clin],ALL.train[sampi,c(16:117)])
	h2=model.matrix(~.,ALL.train2[,-1])
	p.fac2 = rep(1,dim(h2[,-1])[2])
	k=ncol(h2[,-1])-ncol(ALL.train[,c(16:117)])
	p.fac2[c(1:k)] = 0
	XX2<-h2[,-1]
	
	
	res.cv<-try(Fit.Lasso.kCV(X=XX2,y=y,ind.train=ind.train ,ind.test=ind.test,alpha.value=1,penality=p.fac2,indipend.data=indip),silent=T)
	MSE.data[j]<-try(res.cv$MSEi,silent=T)
	MSE.test[j]<-try(res.cv$MCE.test,silent=T)
	MCEtrain[j]<-try(res.cv$MCEtrain,silent=T)
	Pred.test[,j]<-try(as.numeric(res.cv$pred.test),silent=T)
	Pred.train[,j]<-try(as.numeric(res.cv$pred.train),silent=T)
	Pred.data[,j]<-try(as.numeric(res.cv$pred.data),silent=T)
	train[,j]<-ind.train
	test[,j]<-ind.test
	
	print(j)
	
}


save(MCEtrain,MSE.test,MSE.data,train,test,Pred.data,Pred.train,Pred.test,file="3foldCV_LASSO(ALLDATA noPneality)PERMUTATION.Rdata")



nvc<-1000#99 # there was a problem with Cv number 392 
npv.test<-ppv.test<-sen.test<-spe.test<-rep(NA,nvc)
npv.data<-ppv.data<-sen.data<-spe.data<-rep(NA,nvc)
#nvc<-c(1:391,393:1000) # there was a problem with Cv number 392 

for (j in 1:nvc){ 
	
	
	tab.data<-table(y[test[,j]],Pred.data[,j])
	
	npv.data[j]<-tab.data[1,1]/sum(tab.data[1,1],tab.data[2,1])
	ppv.data[j]<-tab.data[2,2]/sum(tab.data[1,2],tab.data[2,2])
	sen.data[j]<-tab.data[2,2]/sum(tab.data[2,1],tab.data[2,2])
	spe.data[j]<-tab.data[1,1]/sum(tab.data[1,2],tab.data[1,1])
	
	tab.test<-table(y2,Pred.test[,j])
	
	npv.test[j]<-tab.test[1,1]/sum(tab.test[1,1],tab.test[2,1])
	ppv.test[j]<-tab.test[2,2]/sum(tab.test[1,2],tab.test[2,2])
	sen.test[j]<-tab.test[2,2]/sum(tab.test[2,1],tab.test[2,2])
	spe.test[j]<-tab.test[1,1]/sum(tab.test[1,2],tab.test[1,1])
	
	print(j)
	
}



save(MCEtrain,MSE.test,MSE.data,train,test,npv.test,ppv.test,sen.test,spe.test,npv.data,
     ppv.data,sen.data,spe.data,file="3foldCV_LASSO(ALLDATA noPneality)PERMUTATION-FINAL.Rdata")


#####################################################

