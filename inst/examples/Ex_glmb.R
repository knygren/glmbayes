set.seed(333)
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
mysd<-1
mu<-matrix(0,5)
mu[1]=log(mean(counts))
V0<-((mysd)^2)*diag(5)
glmb.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu,Sigma=V0))
summary(glmb.D93)

# Menarche Binomial Data Example 
data(menarche2)
Age2=Age-13

#Priors Are Derived in Vignette

## Logit Prior and Model
mu1=as.matrix(c(0,1.098612),ncol=1)
V1<-1*diag(nvars)
V1[1,1]=0.18687882
V1[2,2]=0.10576217
V1[1,2]=-0.03389182
V1[2,1]=-0.03389182
glmb.out1<-glmb(cbind(Menarche, Total-Menarche) ~ Age2,
family=binomial(logit),pfamily=dNormal(mu=mu1,Sigma=V1), data=menarche2)
summary(glmb.out1)

## Probit Prior and Model
mu2=as.matrix(c(0,0.6407758),ncol=1)
V2<-1*diag(nvars)
V2[1,1]=0.07158369
V2[2,2]=0.03453205
V2[1,2]=-0.01512075
V2[2,1]=-0.01512075
glmb.out2<-glmb(cbind(Menarche, Total-Menarche) ~ Age2,
family=binomial(probit),pfamily=dNormal(mu=mu2,Sigma=V2), data=menarche2)
summary(glmb.out2)

## clog-log Prior and Model
mu2=as.matrix(c(-0.3665129,0.6002727),ncol=1)
V2<-1*diag(nvars)
V2[1,1]=0.11491322
V2[2,2]=0.03365986
V2[1,2]=-0.03502783
V2[2,1]=-0.03502783
glmb.out3<-glmb(cbind(Menarche, Total-Menarche) ~ Age2,
                family=binomial(cloglog),pfamily=dNormal(mu=mu2,Sigma=V2), data=menarche2)
summary(glmb.out3)

DIC_Out=rbind(extractAIC(glmb.out1),extractAIC(glmb.out2),extractAIC(glmb.out3))
rownames(DIC_Out)=c("logit","probit","clog-log")
colnames(DIC_Out)=c("pD","DIC")
DIC_Out