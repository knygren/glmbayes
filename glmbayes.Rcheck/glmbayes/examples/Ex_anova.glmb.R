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

# anova for Bayesian Model
# outcome lower Deviance while treament actually increases it
# Use of Traditional anova questionable in Bayesian context
# May Implement Bayes Factor or other approach
anova(glmb.D93)


glm.D93<-glm(counts ~ outcome + treatment, family = poisson())
summary(glm.D93)

# corresponding
anova(glm.D93)