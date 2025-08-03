
## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),pfamily=dNormal(mu=mu,Sigma=V))
## ----glmb confint-------------------------------------------------------------
confint(glmb.D93)
