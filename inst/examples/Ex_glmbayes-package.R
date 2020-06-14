
## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

## ----Prior_and_Calls,results = "hide"-----------------------------------------
## Prior mean vector 
mu<-matrix(0,5)           
mu[1,1]=log(mean(counts)) 
## Prior standard deviation and Variance
mysd<-1           
V=((mysd)^2)*diag(5)  
## Call to glm
glm.D93 <- glm(counts ~ outcome + treatment, 
               family = poisson())
## Call to glmb
glmb.D93<-glmb(n=1000,counts ~ outcome + treatment,
               family = poisson(),mu=mu,Sigma=V)

## ----Printed_Views------------------------------------------------------------
## Printed view of the output from the glm function 
print(glm.D93)
## Printed view of the output from the glmb function 
print(glmb.D93)

## ----glm_Methods--------------------------------------------------------------
## Methods for class "glm""
methods(class="glm")


## ----glmb_Methods-------------------------------------------------------------
## Methods for class "glmb"
methods(class="glmb")

## ----glm_summary--------------------------------------------------------------
## summary output for the "glm" class
summary(glm.D93)

## ----glmb_summary-------------------------------------------------------------
## summary output for the "glm" class
summary(glmb.D93)

## ----glm fitted outputs-------------------------------------------------------
## fitted outputs for the glm function
fitted(glm.D93)

## ----glmb fitted outputs------------------------------------------------------
## mean of fitted outputs for the glm function
## works without a "glmb" class specific generic function
colMeans(fitted(glmb.D93))

## ----glm predictions----------------------------------------------------------
## predictions for the glm function
predict(glm.D93)

## ----glmb predictions---------------------------------------------------------
## predictions for the glmb function
colMeans(glmb.D93$linear.predictors) # no current predict function
colMeans(predict(glmb.D93)) 

## ----glm residuals------------------------------------------------------------
## residuals for the glm function
residuals(glm.D93)

## ----glmb residuals-----------------------------------------------------------
## residuals for the glmb function
colMeans(residuals(glmb.D93))

## ----glm vcov-----------------------------------------------------------------
## nobs for the glm function
vcov(glm.D93)

## ----glmb vcov----------------------------------------------------------------
## nobs for the glm function
vcov(glmb.D93)

## ----glm confint--------------------------------------------------------------
## confint for the glm function
confint(glm.D93)

## ----glmb confint-------------------------------------------------------------
## confint for the glmb function
confint(glmb.D93)


## ----glm AIC------------------------------------------------------------------
## AIC for the glm function (equivalent degrees of freedom and the AIC)
extractAIC(glm.D93)

## ----glmb DIC-----------------------------------------------------------------
## DIC for the glmb function (Estimated effective number of parameters and the DIC)
extractAIC(glmb.D93)

## ----glm Deviance-------------------------------------------------------------
## Deviance for the glm function
deviance(glm.D93)

## ----glmb Deviance------------------------------------------------------------
## Deviance for the glmb function
## works without a "glmb" class specific generic function
mean(deviance(glmb.D93))

## ----glm logLik---------------------------------------------------------------
## Deviance for the glm function
logLik(glm.D93)

## ----glmb logLik--------------------------------------------------------------
## Deviance for the glmb function
mean(logLik(glmb.D93))

## ----glm Model Frame----------------------------------------------------------
## Model Frame for the glm function
model.frame(glm.D93)

## ----glmb Model Frame---------------------------------------------------------
## Model Frame for the glmb function
model.frame(glmb.D93$glm)

## ----glm formula--------------------------------------------------------------
## formula for the glm function
formula(glm.D93)

## ----glmb formula-------------------------------------------------------------
## formula for the glmb function
formula(glmb.D93)

## ----glm nobs-----------------------------------------------------------------
## nobs for the glm function
nobs(glm.D93)

## ----glmb nobs----------------------------------------------------------------
## nobs for the glmb function
nobs(glmb.D93)

## ----glm family---------------------------------------------------------------
## family for the glm function
family(glm.D93)

## ----glmb family--------------------------------------------------------------
## family for the glmb function
family(glmb.D93$glm)

## ----glm show-----------------------------------------------------------------
## nobs for the glm function
show(glm.D93)

## ----glmb show----------------------------------------------------------------
## nobs for the glm function
## works without a "glmb" class specific generic function
show(glmb.D93)


