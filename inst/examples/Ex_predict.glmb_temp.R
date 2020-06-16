
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

## This likely borrow formats from original variables
## 
SF1=matrix(SF[,1],nrow=nrow(newdata1))
SF2=matrix(SF[,2],nrow=nrow(newdata1))
newdata1$SF=cbind(SF1,SF2)
colnames(newdata1$SF)=colnames(SF)
newdata1


pred1=predict(budworm.lgb, newdata=newdata1,olddata=olddata,type = "response",yvars=SF)

predc1=predict(budworm.lg, data.frame(ldose = ld,
sex = factor(rep("M", length(ld)), levels = levels(sex))),type = "response")
predc1
colMeans(pred1)


