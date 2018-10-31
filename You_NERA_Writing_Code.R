# Regression Analysis 
# Huilan You
# The R code of the original paper is included in the first section
# Our modified code is in the second section

# R code for analyzing Bitcoin exchange failure
# Paper http://lyle.smu.edu/~tylerm/fc13.pdf

options(scipen=999)

#read in survival table
survex3<-read.csv("~/Desktop/sdatamodCC.csv", sep=",",header=T)
survex3$dailyvol<-survex3$totalvol/survex3$lifetime

#anti-money laundering indicators
aml<-read.csv("~/Desktop/compliance-aml-cft-whole.csv",sep=",",header=T)
amls<-merge(survex3,aml,by="Country")

amls$PrevN<-amls$Preventive/amls$PreventativeMax
amls$InstN<-amls$Institutional/amls$InstitutionalMax
amls$LegalN<-amls$Legal/amls$LegalMax
amls$AllN<-amls$All/amls$AllMax

amlsv<-amls[!is.na(amls$dailyvol),]
row.names(amlsv)<-amlsv$exchange

library(survival)

# Cox model
cox.vh<-coxph(Surv(time=amlsv$lifetime,event=amlsv$censored,type='right')~log2(amlsv$dailyvol)+amlsv$Hacked+amlsv$All,method="breslow")

summary(cox.vh)

#gives relative risk compared to average
amlsv$risk<-predict(cox.vh,amlsv,type="risk")

plot(survfit(cox.vh))

meancox<-survfit(cox.vh)

meanhaz<-exp(coef(cox.vh)[1]*mean(log2(amlsv$dailyvol),na.rm=T)+coef(cox.vh)[2]*mean(as.numeric(amlsv$Hacked),na.rm=T))

coxplots<-survfit(cox.vh,newdata=amlsv)

amlsv$exchange[!is.na(amlsv$Hacked)&!is.na(amlsv$dailyvol)]

#Code for Figure 1
#pdf("../fig/survcox.pdf",width=6,height=5)
par(mar=c(4.1,4.1,0.5,0.5))
plot(coxplots[15],col="green",lty="dashed",lwd=2,
	xlab="Days",
	ylab="Survival probability",
	cex.lab=1.3,
	cex.axis=1.3						   
)  #Mt Gox					  
lines(coxplots[2],col="red",lty="dotted",lwd=2)		#Vircurex (counterexample: still open,  lo-vol)
lines(coxplots[28],col="blue",lty="dotdash",lwd=2)	#Intersango (still open, higher volume)
lines(coxplots[30],col="brown",lty="longdash",lwd=2)	#Bitfloor (security breach, hard to see what point is).
lines(survfit(cox.vh),lwd=2)				#Mean
legend("topright",legend=c("Intersango","Mt. Gox","Bitfloor","Vircurex","Average"),col=c("blue","green","brown","red","black"),lwd=2,lty=c("dotdash","dashed","longdash","dotted","solid"))
dev.off()

#get coxplots index from the following line (level order)
amlsv$filename


amlsv$lifemo<-amlsv$lifetime/30.#recreate table 1

#Now do the logistic regression based on breach probability
hlogit <- glm(Hacked ~ log2(dailyvol) + lifemo, data = amlsv, family = binomial(link="logit"))
summary(hlogit)

#get the odds ratio with 95% CI
exp(cbind(OR = coef(hlogit), confint(hlogit)))
amlsv$Closed<-amlsv$censor==1

kk=mean(amlsv$lifetime)
#get predictions for hacking probability for hypothetical daily transaction volumes
newdata1 <- with(amlsv, data.frame(lifemo =mean(lifemo), dailyvol=c(10,100,1000,10000)))
newdata1$rankP <- predict(hlogit, newdata = newdata1, type = "response")

#make a graph of predicted hack based on daily volume
newdata2 <- with(amlsv, data.frame(lifemo = mean(lifemo), dailyvol=seq(from = 2, to = 50000, length.out = 10000)))
newdata2 <- cbind(newdata2, predict(hlogit, newdata = newdata2, type = "response", se = TRUE))

newdata2$fitmin<-with(newdata2,fit - (1.645 * se.fit))
newdata2$fitmin[newdata2$fitmin<0]<-0.0
newdata2$fitmax<-with(newdata2,fit + (1.645 * se.fit))
newdata2$fitmax[newdata2$fitmax>1]<-1.0

#Figure 2
#pdf("../fig/volumebreach.pdf",width=6,height=4)
par(mar=c(4.1,4.1,0.5,0.7))
plot(y=newdata2$fit,x=newdata2$dailyvol,type='l',log='x',ylim=c(0,1),lwd=2,
	xlab="Daily transaction volume at exchange",
	ylab="Probability exchange has breach",
	cex.lab=1.3,
	cex.axis=1.3
)
lines(y=newdata2$fitmin,x=newdata2$dailyvol,type='l',lty='dashed',lwd=2)
lines(y=newdata2$fitmax,x=newdata2$dailyvol,type='l',lty='dashed',lwd=2)
legend("topleft",legend=c("Predicted probability","90% C.I."),col="black",lwd=2,lty=c("solid","dashed"))
dev.off()

#stats for model fit
#this gives chi-squared
with(hlogit, null.deviance - deviance)
#this gives associated p value
with(hlogit, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))

# Above is the original code
#------------------------------------------------------------------------
# Below is our own code for 5 different parts of the analysis

#------------------------------------------------------------------------
# 1. alternative model
# A linear regression of the life time of an exchange on daily transaction volume
# the occurence of hacked, and AML compliance
model1=lm(amlsv$lifetime~log2(amlsv$dailyvol)+amlsv$Hacked+amlsv$All,data=amlsv)
summary(model1)

#------------------------------------------------------------------------
# 2. the second alternative model on probability of the closure of the exchange
# A logistic regression model of the probability of the closure on an exchange on 
# daily transaction volume, the occurence of hacked, and AML compliance
logmodel1=glm(amlsv$censored~log2(amlsv$dailyvol)+amlsv$Hacked+amlsv$All,data=amlsv)
summary(logmodel1)

#------------------------------------------------------------------------
# 3. Adding the dummy variable for the North America to the original logistic regression
amlsvNew = amlsv
library("plyr")
amlsvNew <- transform(amlsvNew, Continents=revalue(Continents,c("AS"="REST")))
amlsvNew <- transform(amlsvNew, Continents=revalue(Continents,c("AF"="REST")))
amlsvNew <- transform(amlsvNew, Continents=revalue(Continents,c("AU"="REST")))
amlsvNew <- transform(amlsvNew, Continents=revalue(Continents,c("SA"="REST")))
amlsvNew <- transform(amlsvNew, Continents=revalue(Continents,c("EU"="REST")))
hlogitcontrol <- glm(Hacked ~ log2(dailyvol) + lifemo+Continents, data = amlsvNew, family = binomial(link="logit"))
summary(hlogitcontrol)

#------------------------------------------------------------------------
# 4. Checking the omitted variable 
# dependent variable on independent variable without lifemo
check <- glm(Hacked ~ log2(dailyvol), data = amlsv, family = binomial(link="logit"))
summary(check)
# dependent variable on lifemo
check1 <- glm(Hacked ~ lifemo, data = amlsv, family = binomial(link="logit"))
summary(check1)
# independent variable on lifemo
check2 <- lm(log2(dailyvol)~lifemo, data = amlsv)
summary(check2)
cor(log2(amlsv$dailyvol),amlsv$lifemo)

#------------------------------------------------------------------------
# 5. Test influential points

#Based on the logistic regression of breach probability
#hlogit <- glm(Hacked ~ log2(dailyvol) + lifemo, data = amlsv, family = binomial(link="logit"))
#summary(hlogit)
#Take 1 of the 40 locations out of the list each time and see how the p value on log2(dailyvol) changes

pvlist<-c()
for (i in 1:40){
  Hacked1<-amlsv$Hacked[-i]
  dailyvol1<-amlsv$dailyvol[-i]
  lifemo1<-amlsv$lifemo[-i]
  hlogit1<- glm(Hacked1 ~ log2(dailyvol1) + lifemo1, family = binomial(link="logit"))
  pvlist<-c(pvlist,summary(hlogit1)$coefficients[2,1])
}
print(pvlist)
k<-2
hlogit2<- glm(Hacked[-k] ~ log2(dailyvol[-k]) + lifemo[-k], data = amlsv, family = binomial(link="logit"))
summary(hlogit2)

#------------------------------------------------------------------------






