
require(graphics)

## example from Venables and Ripley (2002, pp. 190-2.)
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
levels(sex)=c("F","M")
SF <- cbind(numdead, numalive = 20-numdead)

budworm.lg <- glm(SF ~ sex*ldose, family = binomial)

################################ Replace this with call to Prior_Set  #########################
## 

#budworm.lg$model
#budworm.lg$prior.weights
betas=budworm.lg$coefficients
mu=matrix(0,nrow=length(betas),ncol=1)
V=sqrt(10)*as.matrix(diag(length(betas)))

rownames(mu)=names(betas)
rownames(V)=names(betas)
colnames(V)=names(betas)
mu
V
##########################################################################

### Run and predict original model
n=1000
budworm.lgb <- glmb(n=n,SF ~ sex*ldose, family = binomial,mu=mu,Sigma=V)
summary(budworm.lgb)
colMeans(predict(budworm.lgb))
olddata=data.frame(SF,sex,ldose)
olddata=data.frame(SF.numdead,SF.numalive,sex,ldose)

##  Add summary function that
#(1) Generates quantiles for predictions
#(2) Simulates for y
#(3) Generates quantiles for y simulations

## For residuals, may want to have summary function that produces quintiles as well

############################################################################

## Generate new data for plot
ld <- seq(0, 5, 0.1)
obs_size=20

##### Male data - Note levels should match those for olddata
newdataM=data.frame(ldose = ld,sex = factor(rep("M", length(ld)), levels = levels(sex)))

## Predict on the response scale (probability)
pred_males=predict(budworm.lgb,newdata=newdataM,olddata=olddata, type="response")

## Simulate new data (also based on 1000 draws from binomial)

pred_males_y=matrix(0,nrow=n,ncol=length(ld))
for(i in 1:n){
  pred_y_males[i,1:length(ld)]=rbinom(length(ld),size=obs_size,prob=pred_males[i,1:length(ld)])
}

## Get means for both

pred_males_m=colMeans(pred_males)
pred_y_males_m=colMeans(pred_y_males/obs_size)

## Get quantiles
quant1_m=apply(pred_males,2,FUN=quantile,probs=c(0.025))
quant2_m=apply(pred_males,2,FUN=quantile,probs=c(0.975))
quant1_m_y=apply(pred_y_males/obs_size,2,FUN=quantile,probs=c(0.025))
quant2_m_y=apply(pred_y_males/obs_size,2,FUN=quantile,probs=c(0.975))

## Set up plot
plot(c(1,32), c(0,1), type = "n", xlab = "dose",
     ylab = "prob", log = "x")

## Plot for males only first 

text(2^ldose, numdead[1:6]/obs_size, label=c("M"))
lines(2^ld, pred_males_m,lty=1)
lines(2^ld, quant1_m,lty=2)
lines(2^ld, quant2_m,lty=2)
lines(2^ld, quant1_m_y,lty=2)
lines(2^ld, quant2_m_y,lty=2)


########################  Female analysis starts here ###########################################


newdataF=data.frame(ldose = ld,sex = factor(rep("F", length(ld)), levels = c("M","F")))
pred_females=predict(budworm.lgb,newdata=newdataF,olddata=olddata, type="response")

## Simulate new data (also based on 20 draws from binomial)

pred_y_females=matrix(0,nrow=n,ncol=length(ld))
for(i in 1:n){
  pred_y_females[i,1:length(ld)]=rbinom(length(ld),size=obs_size,prob=pred_females[i,1:length(ld)])
}

## Get means for both

pred_females_m=colMeans(pred_females)
pred_y_females_m=colMeans(pred_y_females/obs_size)

## Get quantiles
quant1_f=apply(pred_females,2,FUN=quantile,probs=c(0.025))
quant2_f=apply(pred_females,2,FUN=quantile,probs=c(0.975))
quant1_f_y=apply(pred_y_females/obs_size,2,FUN=quantile,probs=c(0.025))
quant2_f_y=apply(pred_y_females/obs_size,2,FUN=quantile,probs=c(0.975))


## Set up plot
plot(c(1,32), c(0,1), type = "n", xlab = "dose",
     ylab = "prob", log = "x")

## Plot for females  

text(2^ldose, numdead[7:12]/obs_size, label=c("F"))
lines(2^ld, pred_females_m,lty=1)
lines(2^ld, quant1_f,lty=2)
lines(2^ld, quant2_f,lty=2)
lines(2^ld, quant1_f_y,lty=2)
lines(2^ld, quant2_f_y,lty=2)
