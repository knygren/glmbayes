## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)

ps=Prior_Setup(weight ~ group)
x=ps$x
mu=ps$mu
V=ps$Sigma
mu[1,1]=mean(weight)

Prior_Check(weight ~ group,family =gaussian(),
            pfamily=dNormal(mu=mu,Sigma=V))
lm.D9 <- lm(weight ~ group,x=TRUE,y=TRUE)
y=lm.D9$y

## Dispersion for maximum likelihood estimate
disp_ML=sigma(lm.D9)^2
n_prior=2
shape=n_prior/2
rate= disp_ML*shape

# Two-Block Gibbs sampler
set.seed(180)
dispersion2=disp_ML
# Run 1000 burn-in iterations
for(i in 1:1000){
  out1=rlmb(n = 1, y=y, x=x, pfamily=dNormal(mu=mu,Sigma=V,dispersion=dispersion2))
  out2=rlmb(n = 1, y=y, x=x, pfamily=dGamma(shape=shape,rate=rate,beta=out1$coefficients[1,]))
  dispersion2=out2$dispersion
}


# Run 1000 additional iterations and store output
beta_out<-matrix(0,nrow=1000, ncol=2)
disp_out=rep(0,1000)
for(i in 1:1000){
  out1=rlmb(n = 1, y=y, x=x, pfamily=dNormal(mu=mu,Sigma=V,dispersion=dispersion2))
  out2=rlmb(n = 1, y=y, x=x, pfamily=dGamma(shape=shape,rate=rate,beta=out1$coefficients[1,]))
  dispersion2=out2$dispersion
  beta_out[i,1:2]=out1$coefficients[1,1:2]
  disp_out[i]=out2$dispersion
}

colMeans(beta_out)
mean(disp_out)

# Same model using Independent_Normal_Gamma_Prior
lmb.D9_v2=lmb(weight ~ group,dIndependent_Normal_Gamma(mu,V,shape=shape,rate=rate))
summary(lmb.D9_v2)
