## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

## Call to glm
glm.D93 <- glm(counts ~ outcome + treatment, 
               family = poisson())

## Using glmb
## Step 1: Set up Prior
ps=Prior_Setup(counts ~ outcome + treatment)
mu=ps$mu
V=ps$Sigma
# Step2A: Check the Prior
Prior_Check(counts ~ outcome + treatment,family = poisson(),
            pfamily=dNormal(mu=mu,Sigma=V))
# Step2B: Update and Re-Check the Prior
mu[1,1]=log(mean(counts))
Prior_Check(counts ~ outcome + treatment,family = poisson(),
            pfamily=dNormal(mu=mu,Sigma=V))
# Step 3: Call the glmb function
glmb.D93<-glmb(counts ~ outcome + treatment, family=poisson(), 
               pfamily=dNormal(mu=mu,Sigma=V))

## ----Printed_Views------------------------------------------------------------
## Printed view of the output from the glm function 
print(glm.D93)
## Printed view of the output from the glmb function 
print(glmb.D93)

## ----Methods---------------------------------------------------------------
## Methods for class "lm"
methods(class="lm")

## Methods for class "glm"
methods(class="glm")

## Methods for class "glmb"
methods(class="glmb")

## ----summary--------------------------------------------------------------
## summary output for the "glm" class
summary(glm.D93)

## summary output for the "glm" class
summary(glmb.D93)

## ----fitted outputs-------------------------------------------------------
## fitted outputs for the glm function
fitted(glm.D93)

## ----glmb fitted outputs------------------------------------------------------
## mean of fitted outputs for the glm function
colMeans(fitted(glmb.D93))

## ----predictions----------------------------------------------------------
## predictions for the glm function
predict(glm.D93)

## predictions for the glmb function
colMeans(predict(glmb.D93)) 

## ----residuals------------------------------------------------------------
## residuals for the glm function
residuals(glm.D93)

## residuals for the glmb function
colMeans(residuals(glmb.D93))

## ----vcov-----------------------------------------------------------------
## vcov for the glm function
vcov(glm.D93)

## vcov for the glm function
vcov(glmb.D93)

## ----confint--------------------------------------------------------------
## confint for the glm function
confint(glm.D93)

## confint for the glm function
confint(glmb.D93)

## ----AIC/DIC------------------------------------------------------------------
## AIC for the glm function (equivalent degrees of freedom and the AIC)
extractAIC(glm.D93)

## DIC for the glmb function (estimated effective number of parameters and the DIC)
extractAIC(glmb.D93)

## ----Deviance-------------------------------------------------------------
## Deviance for the glm function
deviance(glm.D93)

## Deviance for the glmb function
mean(deviance(glmb.D93))

## ----logLik---------------------------------------------------------------
## Deviance for the glm function
logLik(glm.D93)

## Deviance for the glmb function
mean(logLik(glmb.D93))

## ----Model Frame----------------------------------------------------------
## Model Frame for the glm function
model.frame(glm.D93)

## Model Frame for the glmb function
model.frame(glmb.D93$glm)

## ----formula--------------------------------------------------------------
## formula for the glm function
formula(glm.D93)

## ----formula-------------------------------------------------------------
## formula for the glmb function
formula(glmb.D93)

## ----family--------------------------------------------------------------
## family for the glm function
family(glm.D93)

## family for the glmb function
family(glmb.D93$glm)

## ----nobs-----------------------------------------------------------------
## nobs for the glm function
nobs(glm.D93)

## nobs for the glmb function
nobs(glmb.D93)

## ----show-----------------------------------------------------------------
## show for the glm function
show(glm.D93)

## show for the glmb function
show(glmb.D93)
