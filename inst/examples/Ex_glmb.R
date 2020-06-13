## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

n<-1000
mysd<-1

mu<-matrix(0,5)
V0<-((mysd)^2)*diag(5)

glmb.D93<-glmb(n=n,counts ~ outcome + treatment, family = poisson(),mu=mu,Sigma=V0,Gridtype=1)

summary(glmb.D93)

# Menarche Binomial Data Example 

#library(MASS)
data(menarche)

summary(menarche)
plot(Menarche/Total ~ Age, data=menarche)

n<-1000
mu<-matrix(0,nrow=2,ncol=1)
V1<-1*diag(2)
V1[1,1]<-100
V1[2,2]<-1

glm.out1<-glm(cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit), data=menarche)

# temporarily edit out as it seems slow (adding back in - may need
#to adjust prior)

glmb.out1<-glmb(n=n,cbind(Menarche, Total-Menarche) ~ Age,family
                =binomial(logit),mu=mu,Sigma=V1,Gridtype=1, data=menarche)

summary(glm.out1)
summary(glmb.out1)

glm.out2<-glm(cbind(Menarche,Total-Menarche)~Age,family=binomial(probit), data=menarche)
glmb.out2<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age,family=binomial(probit),mu=mu,
                Sigma=V1,Gridtype=1,data=menarche)

summary(glm.out2)
summary(glmb.out2)

glm.out3 <-glm(cbind(Menarche, Total-Menarche) ~ Age,family=binomial(cloglog),
               data=menarche)
glmb.out3<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age,family=binomial(cloglog),
                mu=mu,Sigma=V1,Gridtype=1,data=menarche)

summary(glm.out3)
summary(glmb.out3)