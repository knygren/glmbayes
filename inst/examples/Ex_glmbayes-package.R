## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(),x=TRUE)

n<-1000
mysd<-1

y<-glm.D93$y
x<-glm.D93$x
mu<-matrix(0,5)
#P<-matrix(0.1/((mysd)^2),nrow=5,ncol=5)+1/((mysd)^2)*diag(5)
P<-1/((mysd)^2)*diag(5)

glmtest<-glmb(counts ~ outcome + treatment, family = poisson()
              ,mu=mu,Sigma=solve(P),n=n,Gridtype=1)

glm.D93
glmtest

summary(glm.D93)
summary(glmtest)

fitted(glm.D93)
colMeans(fitted(glmtest))

residuals(glm.D93)
colMeans(residuals(glmtest))

extractAIC(glm.D93)
extractAIC(glmtest)

deviance(glm.D93)
mean(deviance(glmtest))

logLik(glm.D93)
mean(logLik(glmtest))

predict(glm.D93)
#predict(glmtest)  # Need to add function to produce predictions

colMeans(glmtest$linear.predictors)