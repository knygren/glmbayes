data(menarche)

summary(menarche)
plot(Menarche/Total ~ Age, data=menarche)

Age2=menarche$Age-13

x<-matrix(as.numeric(1.0),nrow=length(Age2),ncol=2)
x[,2]=Age2

y=menarche$Menarche/menarche$Total
wt=menarche$Total

mu<-matrix(as.numeric(0.0),nrow=2,ncol=1)
mu[2,1]=(log(0.9/0.1)-log(0.5/0.5))/3

V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates

V1[1,1]<-(0.1*(log(0.9/0.1)-log(0.5/0.5))/2)^2 
V1[2,2]=(0.1*mu[2,1]/2)^2  # Allows slope to be up to 1 times as large as point estimate 

out<-rglmb(n = 1000, y=y, x=x, mu=mu, P=solve(V1), wt = wt, 
family = binomial(logit), Gridtype = 3) 
summary(out)

famfunc<-glmbfamfunc(binomial(logit))

f1<-famfunc$f1
f2<-famfunc$f2
f3<-famfunc$f3
f5<-famfunc$f5
f6<-famfunc$f6

dispersion2<-as.numeric(1.0)
start <- mu
offset2=rep(as.numeric(0.0),length(y))
P=solve(V1)
n=1000

## Appears that the type for some of these arguments are important/problematic

outlist<-glmbsim_NGauss_cpp(n=as.integer(n),y=as.vector(y),
x=as.matrix(x),mu=as.vector(mu),P=as.matrix(P),
offset2=as.vector(offset2),wt=as.vector(wt),dispersion=as.numeric(dispersion2),
famfunc=famfunc,f1=f1,f2=f2,f3=f3,start=as.vector(start),family="binomial",
link="logit",Gridtype=as.integer(3))



### This allows use of the rglmb summary function 
### add interface for glmbsim_NGauss_cpp later

#outlist$call<-match.call()
colnames(outlist$coefficients)<-colnames(x)
class(outlist)<-c(outlist$class,"rglmb")
summary(outlist)

## This envelope seems fine

Env1=outlist$Envelope
Env1

#####

### Adjust weight for dispersion

wt2=wt/dispersion2

################################     #####################

# Shift mean vector to offset so that adjusted model has 0 mean

alpha=x%*%as.vector(mu)+offset2
mu2=0*as.vector(mu)
P2=P
x2=x

############################     ####################################

#####  Optimization step to find posterior mode and associated Precision

parin=start-mu

opt_out=optim(parin,f2,f3,y=as.vector(y),x=as.matrix(x),mu=as.vector(mu),
      P=as.matrix(P),alpha=as.vector(alpha),wt=as.vector(wt2),
      method="BFGS",hessian=TRUE
      )

#opt_out
b2=opt_out$par  ## Posterior mode for adjusted model
b2
b2+as.vector(mu)  # mode for actual model
A1=opt_out$hessian # Approximate Precision at mode


## Try new function here




########################   #################################

# Find eigenvalues and standardize to model with ~ Identity precision at posterior mode

## For consistency with C++ use exported eigenvalue/eigenvector function



A1_eigen=glmb_eig_sym(A1)
D1=diag(A1_eigen$eigval[,1])
L2=sqrt(D1)%*%t(A1_eigen$eigvec)
L2Inv=A1_eigen$eigvec%*%sqrt(solve(D1))



## Apply standardization

b3=L2%*%as.vector(b2)   
#b3  ### Same numbes but ordered differently

b2
A1
A1_eigen$eigval
A1_eigen$eigvec
D1
b3



mu3=L2%*%as.vector(mu2)

x3=x2%*%L2Inv
P3=L2Inv%*%P2%*%t(L2Inv)

P3  # Standardized prior precision - "smaller" than identity matrix

L2Inv%*%b3+as.vector(mu)




#########################            ##########################

###  Find diagonal matrix epsilon so that P3-epsilon is still positive definite

### For this demo, simply use epsilon=0.5*P3
### In C++ code, searches until valid values (multiplying by 0.5 multiple times if needed)

P3Diag=diag(diag(P3))  #   diagonal part of P3
epsilon=0.5*P3Diag
P4=P3-epsilon



epsilon  ## Modified prior precision after P4 shifted to likelihood
P4       ## Precision for Multivariate normal term added to log-likelihood

#########################  Transform model again so that 
###  modified prior is the Standard multivariate normal

A3=diag(length(mu))-epsilon


A3_eigen=glmb_eig_sym(A3)
D2=diag(A3_eigen$eigval[,1])
L3=sqrt(D2)%*%t(A3_eigen$eigvec)
L3Inv=A3_eigen$eigvec%*%sqrt(solve(D2))

A3
D2
A3_eigen


#arma::mat L3= arma::sqrt(D2)*trans(eigvec_2);
#L3Inv=eigvec_2*sqrt(inv_sympd(D2));
b4=L3%*%b3 
b4

mu4=L3%*%mu3 
x4=x3%*%L3Inv
A4=t(L3Inv)%*%A3%*%L3Inv #   Should be transformed data precision matrix
P5=t(L3Inv)%*%P4%*%L3Inv #   Should be precision matrix without epsilon
P6Temp=P5+diag(length(mu)) #  Should be precision matrix for posterior

L3Inv%*%L2Inv%*%b4+as.vector(mu)  ## Check that posterior mode still is the same

#NumericVector bStar, NumericMatrix A, NumericVector y, NumericMatrix x, 
#NumericMatrix mu, NumericMatrix P, NumericVector alpha, NumericVector wt, 
#std::string family, std::string link, int Gridtype, int n, bool sortgrid

Env2=glmbenvelope_c(as.vector(b4), as.matrix(A4),y, as.matrix(x4),
      as.matrix(mu4,ncol=1),as.matrix(P5),as.vector(alpha),as.vector(wt2),
      family="binomial",link="logit",Gridtype=as.integer(3), n=as.integer(n),sortgrid=TRUE)


### Hmmm - this envelope seems wrong  ####

#Envelope
#Env2=outlist$Envelope

Env1$thetabars
Env2$thetabars


