#####################
load("ALLDATA.RData")
ALL.train[,11]<-factor(ALL.train[,11],c("Never","Active","Stopped"))
table(ALL.train[,11])
ALL.test[,11]<-factor(ALL.test[,11],c("Never","Active","Stopped"))
table(ALL.test[,11])
y=ALL.train[,1] ##### THE RESPONSE  TRAIN ####
levels(y)<-c("0","1")
y2=ALL.test[,1] ##### THE RESPONSE  TRAIN ####
levels(y2)<-c("0","1")

NELSON<-ALL.train
NELSON.test<-ALL.test
metab.remove<-c("Var87","Var92","Var96","Var97","Var98","Var99","Var100","Var101","Var102","Var104")
LLP.clin<-NELSON[,c(1:15)]
relevel(LLP.clin[,11],ref=1)#<-c("Never","Active","Stopped ")
LLP.metab<-NELSON[,c(1,16:117)]
LLP.train<-NELSON
LLP.test<-NELSON.test


###### LRT TEST #####
clin<-c("Group","Age","Packyears","BMI","COPD","Cardiac","Smoking","Coagulation") #CHOSEN USING BACKWARD SELECTION AND KEEPING SIGNIFICANT ONES

############## PERMUTING LOGISTIC REGRESSION ####
set.seed(110)
n.cv<-1000#00 # no.of cross validations
MSE.data<-rep(NA,n.cv)
MSE.test<-rep(NA,n.cv)
MCEtrain<-rep(NA,n.cv)
npv.test<-ppv.test<-sen.test<-spe.test<-rep(NA,n.cv)
npv.data<-ppv.data<-sen.data<-spe.data<-rep(NA,n.cv)
for (j in 1:n.cv){ 
	ind.train <-as.vector(sort(sample(1:nrow(LLP.train[,clin]),floor(2*nrow(LLP.train[,clin])/3),replace=F) ) )
	ind.test <- as.vector(c(1:nrow(LLP.train[,clin]))[-ind.train]) 
	logist<-glm(Group~.,data=LLP.train[ind.train,clin],family="binomial")
	pred.fit21 <- try(ifelse(predict(logist,LLP.train[ind.test,clin],type="response")<0.5,"C","LC"),silent=T)
	pred.fit212<- try(ifelse(predict(logist,LLP.test[,-1],type="response")<0.5,"C","LC"),silent=T)
	pred.fit213 <-try(ifelse(predict(logist,LLP.train[ind.train,clin],type="response")<0.5,"C","LC"),silent=T)
	
	MCEtrain[j]<-try(1-sum(diag(table(LLP.train[ind.train,1],pred.fit213)))/length(LLP.train[ind.train,1]),silent=T)
	tab.train<-table(LLP.train[ind.train,1],pred.fit213)
	
	
	MSE.data[j]<-try(1-sum(diag(table(LLP.train[ind.test,1],pred.fit21)))/length(LLP.train[ind.test,1]),silent=T)
	tab.data<-table(LLP.train[ind.test,1],pred.fit21)
	npv.data[j]<-tab.data[1,1]/sum(tab.data[1,1],tab.data[2,1])
	ppv.data[j]<-tab.data[2,2]/sum(tab.data[1,2],tab.data[2,2])
	sen.data[j]<-tab.data[2,2]/sum(tab.data[2,1],tab.data[2,2])
	spe.data[j]<-tab.data[1,1]/sum(tab.data[1,2],tab.data[1,1])
	
	MSE.test[j]<-try(1-sum(diag(table(LLP.test[,1],pred.fit212)))/length(LLP.test[,1]),silent=T)
	tab.test<-table(LLP.test[,1],pred.fit212)
	npv.test[j]<-tab.test[1,1]/sum(tab.test[1,1],tab.test[2,1])
	ppv.test[j]<-tab.test[2,2]/sum(tab.test[1,2],tab.test[2,2])
	sen.test[j]<-tab.test[2,2]/sum(tab.test[2,1],tab.test[2,2])
	spe.test[j]<-tab.test[1,1]/sum(tab.test[1,2],tab.test[1,1])
	
	print(j)
	
}

save(MSE.data,MSE.test,MCEtrain,npv.test,ppv.test,sen.test,spe.test,npv.data,ppv.data,sen.data,spe.data,file="SPLIT-LOGISTIC.RData")
#################
