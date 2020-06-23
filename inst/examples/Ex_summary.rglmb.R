data(menarche2)

summary(menarche2)
plot(Menarche/Total ~ Age, data=menarche2)

Age2=menarche2$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche2$Menarche/menarche2$Total
wt=menarche2$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-((log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(3*mu[2,1]/2)^2  # Allows slope to be up to 3 times as large as point estimate 
prior=list(mu=mu,P=solve(V1))
out<-rglmb(n = 1000, y=y, x=x, prior=prior, wt = wt, 
           family = binomial(logit), Gridtype = 3) 
summary(out)

print(summary(out),digits=4)

mean(out$iters)
