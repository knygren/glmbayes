rm(list=ls())

# Gridtype=1 --> Optimize
# Gridtype=2 --> Use formula to set size?
# Gridtype=3 --> Full size for Grid
# Gridtype=4 --> unidimensional (?)

#Unclear if poor performance because coefficients for two variables essentially zero - Try Alternative
#Problem seems to be primarily if prior means are too far from data --> Leads to prior as "outlier"

## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))


glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(),x=TRUE)

anova(glm.D93)
summary(glm.D93)
predict(glm.D93,newdata=d.AD)


n<-1000
mu<-matrix(0,5)
X<-glm.D93$x
Xmu=glm.D93$x%*%glm.D93$coefficients
explambda=diag(as.vector(exp(Xmu)))

# Use approximately conjugate prior for V0 with 10% weight on prior

wt_0<-0.1
m_0=exp(log(wt_0/(1-wt_0)))
V0<-solve(m_0*t(X)%*%explambda%*%X)

Like_std=summary(glm.D93)$coefficients[,2]
#D93_Prior_Error_Checks=Prior_Likelihood_Check(mu,
#sqrt(diag(V0)),glm.D93$coefficients,Like_std)

#D93_Prior_Error_Checks

#mu[1,1]=0+4*sqrt(diag(V0)[1])
#mu[1,1]=0+15*Like_std[1]
mu[1,1]=log(mean(counts))

#D93_Prior_Error_Checks=Prior_Likelihood_Check(mu,
#sqrt(diag(V0)),glm.D93$coefficients,Like_std)

glmb.D93<-glmb(n=n,counts ~ outcome + treatment, family = poisson(),mu=mu,Sigma=V0,Gridtype=3)




test_dummy=dummy.coef(glmb.D93)

test_dummy[2]

length(test_dummy)


get_all_vars(formula(glm.D93),d.AD)


## This triggers issues with glmb function itsefl

#glmb.D93<-glmb(n=n,counts ~1 , family = poisson(),mu=mu[1],Sigma=V0[1,1],Gridtype=3)

#summary(glmb.D93)

d.AD2 <- data.frame(treatment, outcome)

# This returned a model.frame that looks ok
# First checks are ok

pred_out=predict(glmb.D93,newdata=d.AD2)
colMeans(glmb.D93$linear.predictors)
colMeans(pred_out)

pred_out=predict(glmb.D93,newdata=d.AD2,type="response")
colMeans(glmb.D93$fitted.values)
colMeans(pred_out)


### This likely borrows information from the enviroment
d.AD3 <- data.frame(treatment)
pred_out=predict(glmb.D93,newdata=d.AD3)
colMeans(glmb.D93$linear.predictors)
colMeans(pred_out)


### This should have failed (olddata does not have dependent variable but passed)
### Is rescued dependent variable can be recovered from original frame
pred_out=predict(glmb.D93,newdata=d.AD2)
colMeans(glmb.D93$linear.predictors)
colMeans(pred_out)


### This should fail (wrong number of levels for one of variables in newdata)

outcome2 <- gl(4,1,9)
d.AD4 <- data.frame(treatment, outcome=outcome2, counts)
pred_out=predict(glmb.D93,newdata=d.AD4)  ## Attributes are wrong
pred_out=predict(glmb.D93,newdata=d.AD4,olddata=d.AD) ## Still fails - Attributes are wrong



outcome3 <- gl(3,1,12)
treatment2<-gl(3,3,12)
counts2 <- c(18,17,15,20,10,20,25,13,12,18,17,15)
d.AD5 <- data.frame(treatment=treatment2, outcome=outcome3, counts=counts2)
d.AD6 <- data.frame(treatment=treatment2, outcome=outcome3)

pred_out=predict(glmb.D93,newdata=d.AD5)  # Succeeds because it has counts in it
colMeans(glmb.D93$linear.predictors)
colMeans(pred_out)



pred_out=predict(glmb.D93,newdata=d.AD6)  # Has different lengths from original data - fails
colMeans(glmb.D93$linear.predictors)
colMeans(pred_out)


pred_out=predict(glmb.D93,newdata=d.AD6,olddata=d.AD)  # This now works
colMeans(glmb.D93$linear.predictors)
colMeans(pred_out)


d.AD7=d.AD5[9:12,]
d.AD8=d.AD6[9:12,]

pred_out=predict(glmb.D93,newdata=d.AD7)  # Succeeds because it has counts in it
colMeans(glmb.D93$linear.predictors)
colMeans(pred_out)

pred_out=predict(glmb.D93,newdata=d.AD8)  # Has different lengths from original data - fails
colMeans(glmb.D93$linear.predictors)
colMeans(pred_out)

pred_out=predict(glmb.D93,newdata=d.AD8,olddata=d.AD)  # This works but returns warning
colMeans(glmb.D93$linear.predictors)
colMeans(pred_out)

alpha<-rep(1+rnorm(n=1,mean=0,sd=1),9)

glm.D93_wO<-glm(counts ~ outcome +offset(1+alpha)+ treatment-offset(alpha), family = poisson())
glmb.D93_wO<-glmb(n=n,counts ~ outcome +offset(1+alpha)+ treatment-offset(alpha), family = poisson(),mu=mu,Sigma=V0,Gridtype=3)
summary(glmb.D93_wO)

##### Check out what is goin on here 

pred_out=predict(glmb.D93_wO)

glmb.D93_wO$linear.predictors

colMeans(glmb.D93_wO$linear.predictors)
colMeans(pred_out)

## This one worked
pred_out=predict(glmb.D93_wO,newdata=d.AD2)
colMeans(pred_out)



pred_out=predict(glmb.D93_wO,newdata=d.AD3)
colMeans(pred_out)


pred_out=predict(glmb.D93_wO,newdata=d.AD5)
colMeans(pred_out)


pred_out=predict(glmb.D93_wO,newdata=d.AD5,olddata=d.AD) 
colMeans(pred_out)


d.AD_v2=cbind(d.AD,alpha)
pred_out=predict(glmb.D93_wO,newdata=d.AD5,olddata=d.AD_v2) 
colMeans(pred_out)

d.AD_v2=cbind(d.AD,alpha)
pred_out=predict(glmb.D93_wO,newdata=d.AD5,olddata=d.AD_v2) 
colMeans(pred_out)

d.AD_v2=cbind(d.AD,alpha)
pred_out=predict(glmb.D93_wO,newdata=d.AD7) 
colMeans(pred_out)

### This one should be doable

predict(glm.D93_wO,newdata=d.AD7)

d.AD_v2=cbind(d.AD,alpha)
pred_out=predict(glmb.D93_wO,newdata=d.AD7,olddata=d.AD_v2) 
colMeans(pred_out)



formula(glmb.D93_wO)


get_all_vars(formula(glmb.D93_wO))  ## Variables used to define alpha is in here     

test_form=formula(glmb.D93_wO)
test_terms=terms(test_form)

attr(test_terms,"variables")
attr(test_terms,"variables")[1]
attr(test_terms,"variables")[2]
attr(test_terms,"variables")[3]
attr(test_terms,"variables")[4]
attr(test_terms,"offset") # Columns in model.frame where offsets are located
attr(test_terms,"offset")[1] # offset1
attr(test_terms,"offset")[2] # offset2

length(attr(test_terms,"offset"))

# how to access the offset terms -- need to use terms and the attr offset in order
# to flag which variables in model.frame represent offsets
# the offset terms then need to be summed.
alpha_frame[,attr(test_terms,"offset")[1]]+alpha_frame[,attr(test_terms,"offset")[2]]

## Checking functuon get_all_vars 
## --> Variables and number of rows seem driven by original formula but levels can 
## come from newdata and can contain levels not in original data
## Best to pass without data if goal is to get original variables

get_all_vars(formula(glm.D93))        ## returns 9 observations based on original data 
get_all_vars(formula(glm.D93),d.AD)
get_all_vars(formula(glm.D93),d.AD2)
get_all_vars(formula(glm.D93),d.AD3)
get_all_vars(formula(glm.D93),d.AD4)  ## Contains wrong levels
get_all_vars(formula(glm.D93),d.AD5)  ## Passes
get_all_vars(data=d.AD5)  ## This passes with more rows 

test_frame0=model.frame(formula(glm.D93))
test_frame1=model.frame(formula(glm.D93),d.AD)
test_frame2=model.frame(formula(glm.D93),d.AD2)
test_frame3=model.frame(formula(glm.D93),d.AD3)
test_frame4=model.frame(formula(glm.D93),d.AD4) ## wrong levels included
test_frame5=model.frame(formula(glm.D93),d.AD4,drop.unused.levels = TRUE) ## wrong levels still included
test_frame6=model.frame(formula(glm.D93),d.AD5)  


test_terms2=terms(test_frame2)
attr(test_terms2,"offset") # Columns in model.frame where offsets are located

if(!is.null(attr(test_terms2,"offset"))) print("offset exists")
if(!is.null(attr(test_terms,"offset"))) print("offset exists")


## Fails because lengths differ - Count variable is borrowed from original but had differing length
## from final data
isTRUE(try(model.frame(formula(glm.D93),d.AD6),silent=TRUE))
if(isFALSE(tryCatch(model.frame(formula(glm.D93),d.AD2)))==FALSE) print("This worked")
if(isFALSE(tryCatch(model.frame(formula(glm.D93),d.AD6)))==FALSE) print("This worked")

if(isFALSE(try(model.frame(formula(glm.D93),d.AD2)))==FALSE) print("This worked")
if(isFALSE(try(model.frame(formula(glm.D93),d.AD6)))==FALSE) print("This worked")

testme=try(model.frame(formula(glm.D93),d.AD2))
class(testme)
testme=try(model.frame(formula(glm.D93),d.AD6),silent=TRUE)
class(testme)

isFALSE(tryCatch(model.frame(formula(glm.D93),d.AD6)))

        
        isTRUE(try(model.frame(formula(glm.D93),d.AD6)))

test_frame6=model.frame(formula(glm.D93),d.AD6)  

test_frame6$outcome
attr.all.equal(test_frame1,test_frame2)
attr.all.equal(test_frame1,test_frame4)

test_frame3$outcome
test_frame4$outcome
attr.all.equal(test_frame3$outcome,test_frame4$outcome) # want this to fail
attr.all.equal(test_frame3$outcome,test_frame6$outcome) # want this to pass

names(test_frame4)
nrow(test_frame3)

Compare_Model_Frames(test_frame0,test_frame0)
Compare_Model_Frames(test_frame0,test_frame1)
Compare_Model_Frames(test_frame0,test_frame2)
Compare_Model_Frames(test_frame0,test_frame3)
Compare_Model_Frames(test_frame0,test_frame4)
Compare_Model_Frames(test_frame0,test_frame5)
Compare_Model_Frames(test_frame0,test_frame6)
Compare_Model_Frames(test_frame0,test_frame6,Check_Rows=FALSE)

myerror1=Compare_Model_Frames(test_frame0,test_frame4)
myerror2=Compare_Model_Frames(test_frame0,test_frame6,Check_Rows=FALSE)

model.matrix(formula(glm.D93),test_frame0)

combo_matrix=model.matrix(formula(glm.D93),test_frame6)

## Subsetting loses the attributes
x_old=combo_matrix[1:9,]
x_new=combo_matrix[10:12,]

## Assign composite attributes to sub-matrices
## Maybe this can be applied to the model frame step instead

attr(x_old,which="assign")<-attr(combo_matrix,which="assign")
attr(x_old,which="contrasts")<-attr(combo_matrix,which="contrasts")
x_old

attr(x_new,which="assign")<-attr(combo_matrix,which="assign")
attr(x_new,which="contrasts")<-attr(combo_matrix,which="contrasts")
x_new


## model.frames themselves seem not to have attributes - Only columns do
test_frame0
test_frame6
test_frame6[,names(test_frame6)[1]]
test_frame6[,names(test_frame6)[2]]
test_frame6[,names(test_frame6)[3]]

test_frame7=test_frame6[1:9,]  ## Good news is the attributes for the 
###  levels seem to be retained here -- As long as attributes for composite matches,
##   attributes for the sub-components should match

test_frame7[,names(test_frame7)[1]]
test_frame7[,names(test_frame7)[2]]
test_frame7[,names(test_frame7)[3]]

test_frame8=test_frame6[10:12,]

test_frame8[,names(test_frame8)[1]]
test_frame8[,names(test_frame8)[2]]
test_frame8[,names(test_frame8)[3]]


#attr(test_frame6[,names(test_frame6)[1]],)

attr.all.equal(newframe[, names(newframe)[i]], oldframe[, names(oldframe)[i]])

help(eval)

###################################################################################

head(test_frame4[, c('outcome')])

dim(pred_out$x_old2)
dim(pred_out$x_old)

dim(pred_out$x_old2)[1]

if(!dim(pred_out$x_old2)[1]==dim(pred_out$x_old)[1]) print("Number of rows do not match")
if(!dim(pred_out$x_old2)[2]==dim(pred_out$x_old)[2]) print("Number of columns do not match")


for(i in 1:dim(pred_out$x_old2)[2]){
if(!colnames(pred_out$x_old2)[i]==colnames(pred_out$x_old)[i]) stop("Column names do not match")
}


isTRUE()
isTRUE(all.equal(pred_out$x_old2,pred_out$x_old))

if(isTRUE(all.equal(pred_out$x_old2,pred_out$x_old))==FALSE) print("This is false")


isFALSE(all.equal(pred_out$x_old2,pred_out$x_old))


isTRUE(all.equal(pred_out$x_old,pred_out$x_old))

# Approximate number of candidates per iid sample 1.714
