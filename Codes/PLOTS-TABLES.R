
load("SPLIT-LOGISTIC.RData")

MCE.data.log<-c(paste(round(mean(MSE.data)*100,2),"±",round(sd(MSE.data)*100,2),sep=""))
MCE.test.log<-c(paste(round(mean(MSE.test)*100,2),"±",round(sd(MSE.test)*100,2),sep=""))
SEN.data.log<-c(paste(round(mean(sen.data)*100,2),"±",round(sd(sen.data)*100,2),sep=""))
SEN.test.log<-c(paste(round(mean(sen.test)*100,2),"±",round(sd(sen.test)*100,2),sep=""))
SPE.data.log<-c(paste(round(mean(spe.data)*100,2),"±",round(sd(spe.data)*100,2),sep=""))
SPE.test.log<-c(paste(round(mean(spe.test)*100,2),"±",round(sd(spe.test)*100,2),sep=""))
NPV.data.log<-c(paste(round(mean(npv.data)*100,2),"±",round(sd(npv.data)*100,2),sep=""))
NPV.test.log<-c(paste(round(mean(npv.test)*100,2),"±",round(sd(npv.test)*100,2),sep=""))
PPV.data.log<-c(paste(round(mean(ppv.data)*100,2),"±",round(sd(ppv.data)*100,2),sep=""))
PPV.test.log<-c(paste(round(mean(ppv.test)*100,2),"±",round(sd(ppv.test)*100,2),sep=""))

logistic.data<-c(MCE.data.log,SEN.data.log,SPE.data.log,NPV.data.log,PPV.data.log)
logistic.test<-c(MCE.test.log,SEN.test.log,SPE.test.log,NPV.test.log,PPV.test.log)


load("3foldCV_LASSO(ALLDATA noPneality)FINAL.Rdata")

MCE.data.lasso<-c(paste(round(mean(MSE.data,na.rm=T)*100,2),"±",round(sd(MSE.data,na.rm=T)*100,2),sep=""))
MCE.test.lasso<-c(paste(round(mean(MSE.test,na.rm=T)*100,2),"±",round(sd(MSE.test,na.rm=T)*100,2),sep=""))
SEN.data.lasso<-c(paste(round(mean(sen.data,na.rm=T)*100,2),"±",round(sd(sen.data,na.rm=T)*100,2),sep=""))
SEN.test.lasso<-c(paste(round(mean(sen.test,na.rm=T)*100,2),"±",round(sd(sen.test,na.rm=T)*100,2),sep=""))
SPE.data.lasso<-c(paste(round(mean(spe.data,na.rm=T)*100,2),"±",round(sd(spe.data,na.rm=T)*100,2),sep=""))
SPE.test.lasso<-c(paste(round(mean(spe.test,na.rm=T)*100,2),"±",round(sd(spe.test,na.rm=T)*100,2),sep=""))
NPV.data.lasso<-c(paste(round(mean(npv.data,na.rm=T)*100,2),"±",round(sd(npv.data,na.rm=T)*100,2),sep=""))
NPV.test.lasso<-c(paste(round(mean(npv.test,na.rm=T)*100,2),"±",round(sd(npv.test,na.rm=T)*100,2),sep=""))
PPV.data.lasso<-c(paste(round(mean(ppv.data,na.rm=T)*100,2),"±",round(sd(ppv.data,na.rm=T)*100,2),sep=""))
PPV.test.lasso<-c(paste(round(mean(ppv.test,na.rm=T)*100,2),"±",round(sd(ppv.test,na.rm=T)*100,2),sep=""))

lasso.data<-c(MCE.data.lasso,SEN.data.lasso,SPE.data.lasso,NPV.data.lasso,PPV.data.lasso)
lasso.test<-c(MCE.test.lasso,SEN.test.lasso,SPE.test.lasso,NPV.test.lasso,PPV.test.lasso)

all<-data.frame(rbind(logistic.data,logistic.test,lasso.data,lasso.test))
hed<-c("MCE","SEN","SPE","NPV","PPV")
colnames(all)<-hed
all
library(xtable)
xtable(rbind(hed,logistic.data,logistic.test,lasso.data,lasso.test))

rm(list=ls())
### PLOTS##

par(mfrow=c(1,2))
#### MCE DATA####
load("SPLIT-LOGISTIC.RData")

plot(density(MSE.data,from=0,to=.4,na.rm=T),main="",xlab="MCE",col="blue",ylim=c(0,20),lwd=2)
par(new=T)
load("3foldCV_LASSO(ALLDATA noPneality)FINAL.RData")
plot(density(MSE.data,from=0,to=.4,na.rm=T),main="",xlab="MCE",col="red",ylim=c(0,20),lwd=2)

#### MCE TEST ####
load("SPLIT-LOGISTIC.RData")
plot(density(MSE.test,from=0,to=.4,na.rm=T),main="",xlab="MCE",col="blue",ylim=c(0,30),lwd=2)
par(new=T)
load("3foldCV_LASSO(ALLDATA noPneality)FINAL.RData")
plot(density(MSE.test,from=0,to=.4,na.rm=T),main="",xlab="MCE",col="red",ylim=c(0,30),lwd=2)


######## NPV DATA ####

load("SPLIT-LOGISTIC.RData")
plot(density(npv.data,from=0,to=1,na.rm=T),main=" ",xlab="NPV",col="blue",ylim=c(0,10),lwd=2)
par(new=T)
load("3foldCV_LASSO(ALLDATA noPneality)FINAL.RData")
plot(density(npv.data,from=0,to=1,na.rm=T),main="",xlab="NPV",col="red",ylim=c(0,10),lwd=2)


######## NPV TEST####

load("SPLIT-LOGISTIC.RData")
plot(density(npv.test,from=0,to=1,na.rm=T),main="",xlab="NPV",col="blue",ylim=c(0,20),lwd=2)
par(new=T)
load("3foldCV_LASSO(ALLDATA noPneality)FINAL.RData")
plot(density(npv.test,from=0,to=1,na.rm=T),main=" ",xlab="NPV",col="red",ylim=c(0,20),lwd=2)


#######  PPV DATA#####

load("SPLIT-LOGISTIC.RData")
plot(density(ppv.data,from=0,to=1,na.rm=T),main=" ",xlab="PPV",col="blue",ylim=c(0,10),lwd=2)
par(new=T)
load("3foldCV_LASSO(ALLDATA noPneality)FINAL.RData")
plot(density(ppv.data,from=0,to=1,na.rm=T),main=" ",xlab="PPV",col="red",ylim=c(0,10),lwd=2)

#######  PPV TEST#####

load("SPLIT-LOGISTIC.RData")
plot(density(ppv.test,from=0,to=1,na.rm=T),main=" ",xlab="PPV",col="blue",ylim=c(0,20),lwd=2)
par(new=T)
load("3foldCV_LASSO(ALLDATA noPneality)FINAL.RData")
plot(density(ppv.test,from=0,to=1,na.rm=T),main=" ",xlab="PPV",col="red",ylim=c(0,20),lwd=2)


####### SEN DATA#####  

load("SPLIT-LOGISTIC.RData")
plot(density(sen.data,from=0,to=1,na.rm=T),main=" ",xlab="SENS",col="blue",ylim=c(0,10),lwd=2)
par(new=T)
load("3foldCV_LASSO(ALLDATA noPneality)FINAL.RData")
plot(density(sen.data,from=0,to=1,na.rm=T),main=" ",xlab="SENS",col="red",ylim=c(0,10),lwd=2)

###### SENS TEST#####  

load("SPLIT-LOGISTIC.RData")
plot(density(sen.test,from=0,to=1,na.rm=T),main=" ",xlab="SENS",col="blue",ylim=c(0,15),lwd=2)
par(new=T)
load("3foldCV_LASSO(ALLDATA noPneality)FINAL.RData")
plot(density(sen.test,from=0,to=1,na.rm=T),main=" ",xlab="SENS",col="red",ylim=c(0,15),lwd=2)

##### SPE DATA #####

load("SPLIT-LOGISTIC.RData")
plot(density(spe.data,from=0,to=1,na.rm=T),main=" ",xlab="SPEC",col="blue",ylim=c(0,10),lwd=2)
par(new=T)
load("3foldCV_LASSO(ALLDATA noPneality)FINAL.RData")
plot(density(spe.data,from=0,to=1,na.rm=T),main=" ",xlab="SPEC",col="red",ylim=c(0,10),lwd=2)

##### SPE TEST#####

load("SPLIT-LOGISTIC.RData")
plot(density(spe.test,from=0,to=1,na.rm=T),main=" ",xlab="SPEC",col="blue",ylim=c(0,25),lwd=2)
par(new=T)
load("3foldCV_LASSO(ALLDATA noPneality)FINAL.RData")
plot(density(spe.test,from=0,to=1,na.rm=T),main=" ",xlab="SPEC",col="red",ylim=c(0,25),lwd=2)


######################### PERMUTED DATA #####

load("3foldCV_LASSO(ALLDATA noPneality)PERMUTATION-FINAL.Rdata")
MSE.data1<-as.numeric(MSE.data[-392]) # there was a problem with Cv number 392 
MSE.test1<-as.numeric(MSE.test[-392])
MCE.data.lasso<-c(paste(round(mean(MSE.data1,na.rm=T)*100,2),"±",round(sd(MSE.data1,na.rm=T)*100,2),sep=""))
MCE.test.lasso<-c(paste(round(mean(MSE.test1,na.rm=T)*100,2),"±",round(sd(MSE.test1,na.rm=T)*100,2),sep=""))
SEN.data.lasso<-c(paste(round(mean(sen.data,na.rm=T)*100,2),"±",round(sd(sen.data,na.rm=T)*100,2),sep=""))
SEN.test.lasso<-c(paste(round(mean(sen.test,na.rm=T)*100,2),"±",round(sd(sen.test,na.rm=T)*100,2),sep=""))
SPE.data.lasso<-c(paste(round(mean(spe.data,na.rm=T)*100,2),"±",round(sd(spe.data,na.rm=T)*100,2),sep=""))
SPE.test.lasso<-c(paste(round(mean(spe.test,na.rm=T)*100,2),"±",round(sd(spe.test,na.rm=T)*100,2),sep=""))
NPV.data.lasso<-c(paste(round(mean(npv.data,na.rm=T)*100,2),"±",round(sd(npv.data,na.rm=T)*100,2),sep=""))
NPV.test.lasso<-c(paste(round(mean(npv.test,na.rm=T)*100,2),"±",round(sd(npv.test,na.rm=T)*100,2),sep=""))
PPV.data.lasso<-c(paste(round(mean(ppv.data,na.rm=T)*100,2),"±",round(sd(ppv.data,na.rm=T)*100,2),sep=""))
PPV.test.lasso<-c(paste(round(mean(ppv.test,na.rm=T)*100,2),"±",round(sd(ppv.test,na.rm=T)*100,2),sep=""))

quantile(MSE.data1, probs = c(0.25, 0.975), na.rm =T)
quantile(MSE.test1, probs = c(0.25, 0.975), na.rm =T)

lasso.data<-c(MCE.data.lasso,SEN.data.lasso,SPE.data.lasso,NPV.data.lasso,PPV.data.lasso)
lasso.test<-c(MCE.test.lasso,SEN.test.lasso,SPE.test.lasso,NPV.test.lasso,PPV.test.lasso)

all2<-data.frame(rbind(lasso.data,lasso.test))
hed<-c("MCE","SEN","SPE","NPV","PPV")
colnames(all2)<-hed
all2


