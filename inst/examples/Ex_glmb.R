## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

n<-1000
mysd<-1

mu<-matrix(0,5)
V0<-((mysd)^2)*diag(5)
prior=list(mu=mu,Sigma=V0)
glmb.D93<-glmb(n=n,counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu,Sigma=V0))

summary(glmb.D93)

# Menarche Binomial Data Example 

#library(MASS)
data(menarche2)

summary(menarche2)
plot(Menarche/Total ~ Age, data=menarche2)

n<-1000
mu<-matrix(0,nrow=2,ncol=1)
V1<-1*diag(2)
V1[1,1]<-100
V1[2,2]<-1

glm.out1<-glm(cbind(menarche2$Menarche, menarche2$Total-menarche2$Menarche) ~ menarche2$Age,
              family=binomial(logit), data=menarche2)

# temporarily edit out as it seems slow (adding back in - may need
#to adjust prior)
prior1=list(mu=mu,Sigma=V1)
glmb.out1<-glmb(n=n,cbind(Menarche, Total-Menarche) ~ Age,family
                =binomial(logit),pfamily=dNormal(mu=mu,Sigma=V1), data=menarche2)

summary(glm.out1)
summary(glmb.out1)
glm.out2<-glm(cbind(menarche2$Menarche, menarche2$Total-menarche2$Menarche) ~ menarche2$Age,family=binomial(probit), 
              data=menarche2)
glmb.out2<-glmb(n=n,cbind(menarche2$Menarche, menarche2$Total-menarche2$Menarche) ~ menarche2$Age,
                family=binomial(probit),pfamily=dNormal(mu=mu,Sigma=V1),data=menarche2)

summary(glm.out2)
summary(glmb.out2)

glm.out3 <-glm(cbind(menarche2$Menarche, menarche2$Total-menarche2$Menarche) ~ menarche2$Age,family=binomial(cloglog),
               data=menarche2)
glmb.out3<-glmb(n=n,cbind(menarche2$Menarche, menarche2$Total-menarche2$Menarche) ~ menarche2$Age,family=binomial(cloglog),
                pfamily=dNormal(mu=mu,Sigma=V1))

summary(glm.out3)
summary(glmb.out3)