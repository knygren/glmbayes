library(MASS)
data(menarche)

summary(menarche)
plot(Menarche/Total ~ Age, data=menarche)

n<-1000
mu<-matrix(0,nrow=2,ncol=1)
V1<-1*diag(2)
V1[1,1]<-100
V1[2,2]<-1

glm.out1<-glm(cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit), data=menarche)
glmb.out1<-glmb(n=n,cbind(Menarche, Total-Menarche) ~ Age,family=binomial(logit),mu=mu,Sigma=V1,Gridtype=1, data=menarche)



summary(glm.out1)
summary(glmb.out1)



glm.out2<-glm(cbind(Menarche,Total-Menarche)~Age,family=binomial(probit), data=menarche)
glmb.out2<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age,family=binomial(probit),mu=mu,Sigma=V1,Gridtype=1,data=menarche)

summary(glm.out2)
summary(glmb.out2)




glm.out3 <-glm(cbind(Menarche, Total-Menarche) ~ Age,family=binomial(cloglog), data=menarche)
glmb.out3<-glmb(n=n,cbind(Menarche,Total-Menarche)~Age,family=binomial(cloglog),mu=mu,Sigma=V1,Gridtype=1,data=menarche)



summary(glm.out3)
summary(glmb.out3)




