
# Gridtype=1 --> Optimize
# Gridtype=2 --> Use formula to set size?
# Gridtype=3 --> Full size for Grid
# Gridtype=4 --> unidimensional (?)

#Unclear if poor performance because coefficients for two variables essentially zero - Try Alternative
#Problem seems to be primarily if prior means are too far from data --> Leads to prior as "outlier"

## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))


glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(),x=TRUE)
summary(glm.D93)
predict(glm.D93,newdata=d.AD)


n<-1000
mu<-matrix(0,5)
X<-glm.D93$x
Xmu=glm.D93$x%*%glm.D93$coefficients
explambda=diag(as.vector(exp(Xmu)))

# Use approximately conjugate prior for V0 with 10% weight on prior

wt_0<-0.1
m_0=exp(log(wt_0/(1-wt_0)))
V0<-solve(m_0*t(X)%*%explambda%*%X)

Like_std=summary(glm.D93)$coefficients[,2]
D93_Prior_Error_Checks=Prior_Likelihood_Check(mu,
sqrt(diag(V0)),glm.D93$coefficients,Like_std)

D93_Prior_Error_Checks

#mu[1,1]=0+4*sqrt(diag(V0)[1])
#mu[1,1]=0+15*Like_std[1]
mu[1,1]=log(mean(counts))

D93_Prior_Error_Checks=Prior_Likelihood_Check(mu,
sqrt(diag(V0)),glm.D93$coefficients,Like_std)

glmb.D93<-glmb(n=n,counts ~ outcome + treatment, family = poisson(),mu=mu,Sigma=V0,Gridtype=3)
summary(glmb.D93)

## This triggers issues with glmb function itsefl

#glmb.D93<-glmb(n=n,counts ~1 , family = poisson(),mu=mu[1],Sigma=V0[1,1],Gridtype=3)

#summary(glmb.D93)

d.AD2 <- data.frame(treatment, outcome)
pred_out=predict(glmb.D93,newdata=d.AD2,olddata=d.AD)

colMeans(pred_out)

### This should fail (one of independent variables is missings) and does
d.AD3 <- data.frame(treatment)
pred_out=predict(glmb.D93,newdata=d.AD3,olddata=d.AD)

### This should fail (olddata does not have dependent variable)
pred_out=predict(glmb.D93,newdata=d.AD2,olddata=d.AD2)

### This should fail (wrong number of levels for one of variables in newdata)

outcome2 <- gl(4,1,9)
d.AD4 <- data.frame(treatment, outcome=outcome2, counts)
pred_out=predict(glmb.D93,newdata=d.AD4,olddata=d.AD)



dim(pred_out$x_old2)
dim(pred_out$x_old)

dim(pred_out$x_old2)[1]

if(!dim(pred_out$x_old2)[1]==dim(pred_out$x_old)[1]) print("Number of rows do not match")
if(!dim(pred_out$x_old2)[2]==dim(pred_out$x_old)[2]) print("Number of columns do not match")


for(i in 1:dim(pred_out$x_old2)[2]){
if(!colnames(pred_out$x_old2)[i]==colnames(pred_out$x_old)[i]) stop("Column names do not match")
}


isTRUE()
isTRUE(all.equal(pred_out$x_old2,pred_out$x_old))

if(isTRUE(all.equal(pred_out$x_old2,pred_out$x_old))==FALSE) print("This is false")


isFALSE(all.equal(pred_out$x_old2,pred_out$x_old))


isTRUE(all.equal(pred_out$x_old,pred_out$x_old))

# Approximate number of candidates per iid sample 1.714
