###########################  Example for glmb function:
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

mysd<-1
mu<-matrix(0,5)
mu[1,1]=log(mean(counts))
V0<-((mysd)^2)*diag(5)
glmb.D93<-glmb(counts ~ outcome + treatment, family = poisson(),
pfamily=dNormal(mu=mu,Sigma=V0))
summary(glmb.D93)

###########################  Example for lmb function

## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)

ps=Prior_Setup(weight ~ group)
mu=ps$mu
V=ps$Sigma
mu[1,1]=mean(weight)

Prior_Check(weight ~ group,family =gaussian(),
            pfamily=dNormal(mu=mu,Sigma=V))

## May move this step inside the Prior_Check function
lm.D9 <- lm(weight ~ group,x=TRUE,y=TRUE)
disp_ML=sigma(lm.D9)^2
n_prior=2
shape=n_prior/2
rate= disp_ML*shape

lmb.D9=lmb(weight ~ group,dNormal_Gamma(mu,V/disp_ML,shape=shape,rate=rate))
summary(lmb.D9)
