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

## Will move this step inside the Prior_Check
lm.D9 <- lm(weight ~ group,x=TRUE,y=TRUE)
disp_ML=sigma(lm.D9)^2
n_prior=2
shape=n_prior/2
rate= disp_ML*shape

lmb.D9=lmb(weight ~ group,dNormal_Gamma(mu,V/disp_ML,shape=shape,rate=rate))
summary(lmb.D9)

n<-10000
y<-lmb.D9$y
x<-as.matrix(lmb.D9$x)
prior_list=list(mu=mu,Sigma=V/disp_ML,shape=shape,rate=rate)
ngamma.D9=rNormal_Gamma_reg(n=1000,y=y,x=x,
  prior_list=prior_list)
summary(ngamma.D9$coefficients)

