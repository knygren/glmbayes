## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))


glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(),x=TRUE)

summary(glm.D93)

n<-1000
mu<-matrix(0,5)
X<-glm.D93$x
Xmu=glm.D93$x%*%glm.D93$coefficients
explambda=diag(as.vector(exp(Xmu)))

# Use approximately conjugate prior for mu with 10% weight on prior

wt_0<-0.1
m_0=exp(log(wt_0/(1-wt_0)))
V0<-solve(m_0*t(X)%*%explambda%*%X)

diag(sqrt(diag(V0)))
V0
sqrt(diag(V0))


Like_std=summary(glm.D93)$coefficients[,2]
D93_Prior_Error_Checks=Prior_Likelihood_Check(glm.D93,mu,V0)

D93_Prior_Error_Checks


mu[1,1]=0+4*sqrt(diag(V0)[1])
mu[1,1]=0+15*Like_std[1]
mu[1,1]=log(mean(counts))
mu

D93_Prior_Error_Checks=Prior_Likelihood_Check(glm.D93,mu,V0)

#glmb.D93<-glmb(n=n,counts ~ outcome + treatment, family = poisson(),mu=mu,Sigma=V0,Gridtype=3)
#summary(glmb.D93)
