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
glmb.D93<-glmb(n=n,counts ~ outcome + treatment, family = poisson(),prior=prior)

summary(glmb.D93)