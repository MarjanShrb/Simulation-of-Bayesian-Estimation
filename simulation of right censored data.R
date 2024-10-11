#-------------------------------------------------------------------------------
###---------------------------simulation of Bayesian estimation of 
###---------------------------proportional hazards model parameters 
###---------------------------with length-biased and right censored data--------
#-------------------------------------------------------------------------------
rm(list = ls())

#--------beta=c(1,1),light censoring---------------

library(Rlab)
n=400
set.seed(13761376)
Z=matrix(0,n,2)
T=numeric(n)
A=numeric(n)
k=1
set.seed(13761376)
while(k<(n+1)) {
  x1=rbern(1,0.5)
  x2=runif(1,-0.5,0.5)
  t=sqrt((-log(runif(1)))/(exp(x1+0.2*x2)))
  a=runif(1,0,2.2)
  if(t>a){
    Z[k,1]=x1
    Z[k,2]=x2
    T[k]=t
    A[k]=a
    k=k+1
  }
}
#censoring 
set.seed(13761376)
C=runif(n,0,0.5)
delta=numeric(n)
Y=numeric(n)
set.seed(13761376)
for (i in 1:n) {
  if(T[i]>C[i]+A[i]){
    Y[i]=A[i]+C[i]
    delta[i]=0
  }
  else {
    Y[i]=T[i]
    delta[i]=1
  }
}
.

# knots and gamma ---------------------------------------------------------

knots_A=c(0)
for (i in 1:(5+1)) {
  knots_A[i+1]=i*((2)/(6))
}
knots_Y=c(0)
for (i in 1:(5+1)) {
  knots_Y[i+1]=i*((2)/(6))
}
m=7
k=2
s_Y=c(rep(1,k)*knots_Y[1], knots_Y[2:(m-1)], rep(1,k)*knots_Y[m])
s_A=c(rep(1,k)*knots_A[1], knots_A[2:(m-1)], rep(1,k)*knots_A[m])
s=s_A

###--------------------------posterior d=1 full likelihood
library(splines2)
library(Rmpfr)
pos_beta1_d1=function(be1){
  final=c()
  IsY=c()
  for(i in 1:n){ 
    IsY[i]=sum(gamma*iSpline(Y[i],knots =knots_Y[-c(1,7)],degree = 1,Boundary.knots=c(0,2)))
  }
  for (i in 1:n) {
    DEN=c(0)
    for (l in 2:7) {
      lp=sum(gamma[1:(l-1)])
      DEN[l]=exp( -exp(be1*Z[i,1]+be2*Z[i,2])*lp )*
        integrate( Vectorize(function(a) exp( ( exp(be1*Z[i,1]+be2*Z[i,2])*(gamma[l-1]*(s[l+1]-a)^2 *(s[l+2]-s[l])-
                                                                  gamma[l]*(a-s[l])^2 *(s[l+1]-s[l-1])) )/
                                                ( (s[l+1]-s[l])*(s[l+1]-s[l-1])*(s[l+2]-s[l]) ) )) ,s[l],s[l+1] )$value
    }
    
    final[i]= (exp(-IsY*exp(be1*Z[i,1]+be2*Z[i,2]))*exp((be1*Z[i,1]+be2*Z[i,2])*delta[i]))/
      (sum(DEN))
  }
  
  prod(mpfr(final,2))*dnorm(be1,0,10)
}


pos_beta2_d1=function(be2){
  final=c()
  IsY=c()
  for(i in 1:n){ 
    IsY[i]=sum(gamma*iSpline(Y[i],knots =knots_Y[-c(1,7)],degree = 1,Boundary.knots=c(0,2)))
  }
  for (i in 1:n) {
    DEN=c(0)
    for (l in 2:7) {
      lp=sum(gamma[1:(l-1)])
      DEN[l]=exp( -exp(be1*Z[i,1]+be2*Z[i,2])*lp )*
        integrate( Vectorize(function(a) exp( ( exp(be1*Z[i,1]+be2*Z[i,2])*(gamma[l-1]*(s[l+1]-a)^2 *(s[l+2]-s[l])-
                                                                              gamma[l]*(a-s[l])^2 *(s[l+1]-s[l-1])) )/
                                                ( (s[l+1]-s[l])*(s[l+1]-s[l-1])*(s[l+2]-s[l]) ) )) ,s[l],s[l+1] )$value
    }
    
    final[i]= (exp(-IsY*exp(be1*Z[i,1]+be2*Z[i,2]))*exp((be1*Z[i,1]+be2*Z[i,2])*delta[i]))/
      (sum(DEN))
  }
  
  prod(mpfr(final,2))*dnorm(be2,0,10)
}


# requirements ------------------------------------------------------------
set.seed(13761376)
#lambda=rgamma(1,1,1)
lambda=1
set.seed(13761376)
gamma=rexp(7,lambda)


##------------------ARS package------------------------------

library(armspp)

b1pos=function(x)as.numeric(log(pos_beta1_d1(x)))
b2pos=function(x)as.numeric(log(pos_beta2_d1(x)))

be2=0.2
mean(arms(100,b1pos,0.9,4,metropolis = TRUE))
be1=1
mean(arms(100,b2pos,-4,4,metropolis = TRUE))






iteration=10

beta1=matrix(0,iteration,100)
simbeta1=c()
rm(be1)
be2=0.2
for(i in 1:iteration){
  set.seed(as.integer(Sys.time()))
  beta1[i,]=arms(50,b1pos,0.9,4,metropolis = TRUE)
  simbeta1[i]=mean(beta1[i,])
}
beep(8)
bias1=round(1-mean(simbeta1),4)


simbeta2=c()
beta2=matrix(0,iteration,100)
rm(be2)
be1=1
for(i in 1:iteration){
  set.seed(as.integer(Sys.time()))
  beta2[i,]=arms(50,b2pos,-4,4,metropolis = TRUE)
  simbeta2[i]=mean(beta2[i,])
}
beep(8)
bias2=round(0.2-mean(simbeta2),4)

SSD1=round(sqrt((sum((simbeta1-mean(simbeta1))^2))/(9)),4)
SSD2=round(sqrt((sum((simbeta2-mean(simbeta2))^2))/(9)),4)

ESMSE1=round(sqrt(bias1^2+SSD1^2),4)
ESMSE2=round(sqrt(bias2^2+SSD2^2),4)

in_interval=c()
for (i in 1:iteration) {
  ci_lower <- quantile(beta1[i,], 0.035)
  ci_upper <- quantile(beta1[i,], 0.855)
  in_interval[i] <- (1 >= ci_lower) && (1 <= ci_upper)
}
mean(in_interval)

in_interval=c()
for (i in 1:iteration) {
  ci_lower <- quantile(beta2[i,], 0.025)
  ci_upper <- quantile(beta2[i,], 0.975)
  in_interval[i] <- (0.2 >= ci_lower) && (0.2 <= ci_upper)
}
mean(in_interval)


##-------------------cpxph-----------

library(CoxPhLb)
aa=data.frame(A,Y,delta,Z)
fit= coxphlb(Surv(A,Y,delta)~ (X1 + X2),data=aa,method = "EE")
bias1=1-coef(fit)[1]
bias2=0.2-coef(fit)[2]

SSD1=fit$std.err[1]
SSD2=fit$std.err[2]

ESMSE1=round(sqrt(bias1^2+SSD1^2),4)
ESMSE2=round(sqrt(bias2^2+SSD2^2),4)

in_interval=c()
for (i in 1:iteration) {
  ci_lower <- quantile(beta1[i,], 0.035)
  ci_upper <- quantile(beta1[i,], 0.855)
  in_interval[i] <- (1 >= ci_lower) && (1 <= ci_upper)
}
mean(in_interval)

in_interval=c()
for (i in 1:iteration) {
  ci_lower <- quantile(beta2[i,], 0.025)
  ci_upper <- quantile(beta2[i,], 0.975)
  in_interval[i] <- (0.2 >= ci_lower) && (0.2 <= ci_upper)
}
mean(in_interval)

time_frame <- 100
predicted_survival_probabilities <- predict(fit,newdata  =aa,type  ="survival")
coverage_probability<- sum( data$event  == 1&data$ time  <= time_frame)/ nrow( data)
