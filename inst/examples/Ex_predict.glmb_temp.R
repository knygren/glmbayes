
require(graphics)

## example from Venables and Ripley (2002, pp. 190-2.)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)

ldose <- rep(0:5, 2)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive = 20-numdead)

olddata=data.frame(ldose = ldose,sex = factor(sex,levels=c("F","M")))
olddata$SF=SF
olddata  ## Names here are correct

## Add SF here instead

#olddata=data.frame(ldose = ldose,sex = factor(sex,levels=c("F","M")))

budworm.lg <- glm(SF ~ sex*ldose, family = binomial,x=TRUE,data=olddata)
budworm.lgb <- glmb(n=1000,SF ~ sex*ldose, family = binomial,mu=matrix(0,nrow=4,ncol=1),
                    Sigma=10*diag(4),data=olddata)

ld <- seq(0, 5, 0.1)
newdata1=data.frame(ldose = ld, sex = factor(rep("M", length(ld)), levels = c("F","M")))

yvars=as.character(formula(budworm.lgb))[2]
oldyvars=olddata[yvars]
m_oldyvars=colMeans(oldyvars)
names(m_oldyvars)[1]

names(budworm.lgb$glm$model)
colnames(oldyvars[,1])

rm(SF_Temp2)

for(i in 1:length(names(m_oldyvars)))
  
{
Temp1=matrix(colMeans(olddata[yvars])[i],nrow=nrow(newdata1),ncol=1)  
colnames(Temp1)=names(m_oldyvars)[i]
if(i==1) Temp2=Temp1
else(Temp2=cbind(Temp2,Temp1))
  
}

Temp2

SF_Temp2_1=matrix(colMeans(olddata[yvars])[1],nrow=nrow(newdata1),ncol=1)
colnames(SF_Temp2_1)=names(m_oldyvars)[1]
SF_Temp2_2=matrix(colMeans(olddata[yvars])[2],nrow=nrow(newdata1),ncol=1)
colnames(SF_Temp2_2)=names(m_oldyvars)[2]

SF_Temp=SF_Temp2_1
SF_Temp=cbind(SF_Temp,SF_Temp2_2)
SF_Temp


SF_Temp[1]=rep(colMeans(olddata[yvars])[1],nrow(newdata1))




newdata1$SF=matrix(colMeans(olddata[yvars]),nrow=nrow(newdata1),
                   ncol=length(colMeans(olddata[yvars])),byrow = TRUE)
colnames(newdata1$SF)=colnames(oldyvars[,1])
olddata




## This likely borrow formats from original variables
## 
SF1=matrix(SF[,1],nrow=nrow(newdata1))
SF2=matrix(SF[,2],nrow=nrow(newdata1))
newdata1$SF=cbind(SF1,SF2)
colnames(newdata1$SF)=colnames(SF)
newdata1



pred1=predict(budworm.lgb, newdata=newdata1,olddata=olddata,type = "response",yvars=SF)

olddatatemp=olddata[c("ldose","sex")]
pred1=predict(budworm.lgb, newdata=newdata1,olddata=olddatatemp,type = "response",yvars=SF)

newdatatemp=newdata1[c("SF","sex")]
pred1=predict(budworm.lgb, newdata=newdatatemp,olddata=olddata,type = "response",yvars=SF)

newdatatemp=newdata1[c("ldose","sex")]
pred1=predict(budworm.lgb, newdata=newdatatemp,olddata=olddata,type = "response",yvars=SF)

colMeans(pred1)


predc1=predict(budworm.lg, data.frame(ldose = ld,
sex = factor(rep("M", length(ld)), levels = levels(sex))),type = "response")
predc1
colMeans(pred1)

## Use this to extract all variables from original model formula
mod_vars=all.vars(formula(budworm.lgb))

yvars=as.character(formula(budworm.lgb))[2]
xvars <- setdiff(mod_vars, yvars)
xvars
yvars

temp_vars=c(mod_vars,"no_good")
olddata[mod_vars]

tryCatch(olddata[temp_vars],error=function(e) {stop("olddata does not contain all model variables")})


rm(list=ls())

as.character(formula(budworm.lgb))[1]
yvars=as.character(formula(budworm.lgb))[2]
as.character(formula(budworm.lgb))[3]

olddata[yvars]

names(olddata[,1])
colMeans(olddata[yvars])



## 


