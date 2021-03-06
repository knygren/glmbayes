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

## Start setup here [First output from earlier optim optimization]
betastar=glmb.D93$coef.mode  # Posterior mode from optim
x=glmb.D93$x
y=glmb.D93$y
#offset=glmb.D93$offset   # not present in the output --> For now set to 0 vector
offset2=0*y   # Should return this from lower level functions
weights2=glmb.D93$prior.weights

## Check influence measures for original model
fit=glmb.wfit(x,y,weights2,offset2,family=poisson(),Bbar=mu,P=solve(V0),betastar)
influence.measures(fit)

print(fit)
print(glmb.D93$coef.mode)

### Now try a strong prior with poorly chosen intercept
mu1=0*mu
V1=0.1*V0
glmb2.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=mu1,Sigma=V1))

Bbar2=mu1  # Prior mean
betastar2=glmb2.D93$coef.mode  # Posterior mode from optim
fit2=glmb.wfit(x,y,weights2,offset2,family=poisson(),Bbar2,P=solve(V1),betastar2)

influence.measures(fit2)

print(fit2)
print(glmb2.D93$coef.mode)



influence(glmb2.D93)
influence(glmb2.D93,do.coef=TRUE)
influence(glmb2.D93,do.coef=FALSE)

#glmb2.D93$qr=glmb2.D93$fit$qr
#influence.measures(glmb2.D93,influence(glmb2.D93))
#influence.measures(glmb2.D93$fit,influence(glmb2.D93))
#influence.measures(glmb2.D93$fit,influence(glmb2.D93,do.coef=FALSE))
glmb.influence.measures(glmb2.D93,influence(glmb2.D93))


## Do not need method functions
dfbeta(glmb2.D93,influence(glmb2.D93))
dfbeta(glmb2.D93$fit,influence(glmb2.D93))
hatvalues(glmb2.D93,influence(glmb2.D93))
hatvalues(glmb2.D93$fit,influence(glmb2.D93))

# Methods (implmented)

rstandard(glmb2.D93,infl=influence(glmb2.D93))
rstandard(glmb2.D93)


# Need Metods

## Needs a method function
#dfbetas(glmb2.D93,influence(glmb2.D93))
dfbetas(glmb2.D93)
dfbetas(glmb2.D93,influence(glmb2.D93,do.coef=TRUE))
dfbetas(glmb2.D93$fit,influence(glmb2.D93))




# Needs a method function
cooks.distance(glmb2.D93)
cooks.distance(glmb2.D93$fit,influence(glmb2.D93))



# Needs a method function (now works)
#rstandard(glmb2.D93$fit,influence(glmb2.D93))


# Needs a method function
rstudent(glmb2.D93)
rstudent(glmb2.D93$fit,influence(glmb2.D93))


# Not a method, separate function
glmb.dffits(glmb2.D93)
dffits(glmb2.D93$fit,influence(glmb2.D93)) 

# Not a method - separate function
glmb.covratio(glmb2.D93)  
covratio(glmb2.D93$fit,influence(glmb2.D93))  



# Needs a method function but requires a different approach
# The influence function actually stores this measure
## This seems to require more work to create a method function
## Should 
#hat(glmb2.D93,influence(glmb2.D93))
#hat(glmb2.D93$fit,influence(glmb2.D93))

infl=influence(glmb2.D93)
hat2=infl$hat               
hat2               
