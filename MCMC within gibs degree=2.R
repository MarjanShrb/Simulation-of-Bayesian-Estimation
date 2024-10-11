##------------------------------------------------------------------------------
##----------------MCMC within Gibbs method by parameter(b1,b2)(gamma) ----------
###-------------------------simulation d=2  without integral -------------------
##-------------------------- sample size----------------------------------------
##------------------------------------------------------------------------------

##------------------------data generate----------------------------
rm(list = ls())
library(Rlab)
n=25
set.seed(13761376)
Z=matrix(0,n,2)
U <- vector(mode = "list", length = n)
T=numeric(n)
A=numeric(n)
LR=matrix(0,n,4)
k=1
set.seed(13761376)
while(k<(n+1)) {
  z1=rbern(1,0.5)
  z2=runif(1,-0.5,0.5)
  t=(log(runif(1)))/(-0.3*exp(z1+z2))
  a=runif(1,0,14)
  if(t>=a){
    Z[k,1]=z1
    Z[k,2]=z2
    T[k]=t
    A[k]=a
    U[[k]][1]=A[k]
    i=1
    while(U[[k]][i]<=14){
      uu=U[[k]][i]+0.1+runif(1,0,2)
      if(uu<=14) {
        U[[k]][i+1]=uu
        i=i+1
      }
      else{
        break
      }
    }
    i=i-1
    for(j in 1:i){
      if(U[[k]][j]<T[k] & T[k] <= U[[k]][j+1]) {
        LR[k,1]=U[[k]][j]
        LR[k,2]=U[[k]][j+1]
        LR[k,3]=j
        LR[k,4]=j+1
      }
      else if(T[k]>U[[k]][i+1]){
        LR[k,1]=U[[k]][i+1]
        LR[k,2]=Inf
        LR[k,3]=i+1
        LR[k,4]=Inf
      }
    }
    k=k+1
  }
}


#----------------------posterior def with integral ---------------------

library(splines2)

post_beta1_d1=Vectorize(function(be1){
  IsL=c()
  IsR=c()
  for(i in 1:n){ 
    IsL[i]=sum(gamma*iSpline(LR[i,1], knots = knots[-c(1,22)], degree = 1, Boundary.knots=c(0,14)))
    IsR[i]=sum(gamma*iSpline(LR[i,2], knots = knots[-c(1,22)], degree = 1, Boundary.knots=c(0,14)))
  }
  final=c()
  for(i in 1:n){
    DEN=c(0)
    for (l in 2:(interior_knot)) {
      lp=sum(gamma[1:(l-1)])
      DEN[l]=exp( -exp(be1*Z[i,1]+be2*Z[i,2])*lp )*
        integrate( Vectorize(function(a) exp( ( exp(be1*Z[i,1]+be2*Z[i,2])*(gamma[l-1]*(s[l+1]-a)^2 *(s[l+2]-s[l])-
                                                                              gamma[l]*(a-s[l])^2 *(s[l+1]-s[l-1])) )/
                                                ( (s[l+1]-s[l])*(s[l+1]-s[l-1])*(s[l+2]-s[l]) ) )) ,s[l],s[l+1] )$value
    }
    final[i]=(exp(-(IsL[i])*exp(be1*Z[i,1]+be2*Z[i,2]))-
                (ifelse(LR[i,2]<Inf,exp(-( IsR[i])*exp(be1*Z[i,1]+be2*Z[i,2])),0)))/
      ((sum(DEN))/14)
  }
  prod(final)*dnorm(be1,1,5)
})

post_beta2_d1=Vectorize(function(be2){
  IsL=c()
  IsR=c()
  for(i in 1:n){ 
    IsL[i]=sum(gamma*iSpline(LR[i,1], knots = knots[-c(1,22)], degree = 1, Boundary.knots=c(0,14)))
    IsR[i]=sum(gamma*iSpline(LR[i,2], knots = knots[-c(1,22)], degree = 1, Boundary.knots=c(0,14)))
  }
  final=c()
  for(i in 1:n){
    DEN=c(0)
    for (l in 2:(interior_knot)) {
      lp=sum(gamma[1:(l-1)])
      DEN[l]=exp( -exp(be1*Z[i,1]+be2*Z[i,2])*lp )*
        integrate( Vectorize(function(a) exp( ( exp(be1*Z[i,1]+be2*Z[i,2])*(gamma[l-1]*(s[l+1]-a)^2 *(s[l+2]-s[l])-
                                                                              gamma[l]*(a-s[l])^2 *(s[l+1]-s[l-1])) )/
                                                ( (s[l+1]-s[l])*(s[l+1]-s[l-1])*(s[l+2]-s[l]) ) )) ,s[l],s[l+1] )$value
    }
    final[i]=(exp(-(IsL[i])*exp(be1*Z[i,1]+be2*Z[i,2]))-
                (ifelse(LR[i,2]<Inf,exp(-( IsR[i])*exp(be1*Z[i,1]+be2*Z[i,2])),0)))/
      ((sum(DEN))/14)
  }
  prod(final)*dnorm(be2,1,5)
})

post_gamma_d1=Vectorize(function(dxp,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,
                                 g18,g19,g20,g21){
  Isgamma_L=c()
  Isgamma_R=c()
  gamma=c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,
          g18,g19,g20,g21)
  for(i in 1:n){ 
    Isgamma_L[i]=sum(gamma*iSpline(LR[i,1], knots = knots[-c(1,22)], degree = 1, Boundary.knots=c(0,14)))
    Isgamma_R[i]=sum(gamma*iSpline(LR[i,2], knots = knots[-c(1,22)], degree = 1, Boundary.knots=c(0,14)))
  }
  final=c()
  for(i in 1:n){
    DEN=c(0)
    for (l in 2:(interior_knot)) {
      lp=sum(gamma[1:(l-1)])
      DEN[l]=exp( -exp(be1*Z[i,1]+be2*Z[i,2])*lp )*
        integrate( Vectorize(function(a) exp( ( exp(be1*Z[i,1]+be2*Z[i,2])*(gamma[l-1]*(s[l+1]-a)^2 *(s[l+2]-s[l])-
                                                                              gamma[l]*(a-s[l])^2 *(s[l+1]-s[l-1])) )/
                                                ( (s[l+1]-s[l])*(s[l+1]-s[l-1])*(s[l+2]-s[l]) ) )) ,s[l],s[l+1] )$value
    }
    final[i]=(exp(-(Isgamma_L[i])*exp(be1*Z[i,1]+be2*Z[i,2]))-
                (ifelse(LR[i,2]<Inf,exp(-(Isgamma_R[i])*exp(be1*Z[i,1]+be2*Z[i,2])),0)))/
      ((sum(DEN))/14)
  }
  prod(final)*dexp(dxp,lamda_hat)
})



#------------------------------simulation ----------------------------------

##---------requirments
set.seed(13761376)
lamda_hat=rgamma(1,1,1)
interior_knot=20
gamma=numeric(interior_knot+1)
set.seed(13761376)
for (i in 1:(interior_knot+1)) {
  gamma[i]=rexp(1,lamda_hat)
}
knots=c(0)
for (i in 1:(interior_knot+1)) {
  knots[i+1]=i*((14)/(21))
}
#knots=round(knots,2)
m=length(knots)
k=1
s=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m])


kk=c(0,0)
be2=0
xbeta1=c()
xbeta1[1] <- runif(1,0.5,1.5)
xbeta2=c()
xbeta2[1] <- runif(1,0.5,1.5)

####----------------------------beta 1 and beta2
for (i in 2:10000) {
  
  ###----------------------------Beta 1
  u <- runif(1)
  xt <- xbeta1[i-1]
  y <- runif(1,xt-1,xt+1)
  num <- post_beta1_d1(y) * dunif(xt, y-1,y+1)
  den <- post_beta1_d1(xt) * dunif(y, xt-1,xt+1)
  if (u <= num/den) xbeta1[i] <- y 
  else {
    xbeta1[i] <- xt
    kk[1] <- kk[1]+1 #y is rejected
  }
  
  ###-----------------------------Beta 2
  be1=xbeta1[i]
  u <- runif(1)
  xt <- xbeta2[i-1]
  y <- runif(1,xt-1,xt+1)
  num <- post_beta2_d1(y) * dunif(xt, y-1,y+1)
  den <- post_beta2_d1(xt) * dunif(y, xt-1,xt+1)
  if (u <= num/den) xbeta2[i] <- y 
  else {
    xbeta2[i] <- xt
    kk[2] <- kk[2]+1 #y is rejected
  }
  
  ###next step  
  be2=xbeta2[i]
  
}


###gamma----------------------------
be1=mean(xbeta1[1000:10000])
be2=mean(xbeta2[1000:10000])
xgamma=matrix(0,10000,21)
for (i in 1:21) {
  xgamma[1,i]=rexp(1,lamda_hat)
}
kg=rep(0,21)
vog= vector(mode = "list", length = 21)
for (i in 1:21) {
  vog[[i]]=gamma[i]
}
vog_x=vector(mode = "list", length = 21)
for (i in 1:21) {
  vog_x[[i]]=gamma[i]
}


for (i in 2:10000) {
  for (j in 1:21) {
    u <- runif(1)
    xt <- xgamma[i-1,j]
    y <- rexp(1,xt+1)
    vog[[j]]=y
    vog_x[[j]]=xt
    posgy=post_gamma_d1(vog[[j]],vog[[1]],vog[[2]],vog[[3]],vog[[4]],vog[[5]],vog[[6]],
                        vog[[7]],vog[[8]],vog[[9]],vog[[10]],vog[[11]],vog[[12]],
                        vog[[13]],vog[[14]],vog[[15]],vog[[16]],vog[[17]],vog[[18]],
                        vog[[19]],vog[[20]],vog[[21]] )
    posgxt=post_gamma_d1(vog_x[[j]],vog_x[[1]],vog_x[[2]],vog_x[[3]],vog_x[[4]],vog_x[[5]],vog_x[[6]],
                         vog_x[[7]],vog_x[[8]],vog_x[[9]],vog_x[[10]],vog_x[[11]],vog_x[[12]],
                         vog_x[[13]],vog_x[[14]],vog_x[[15]],vog_x[[16]],vog_x[[17]],vog_x[[18]],
                         vog_x[[19]],vog_x[[20]],vog_x[[21]] )
    num <- posgy * dexp(xt,y)
    den <- posgxt * dexp(y,xt)
    if (u <= num/den) xgamma[i,j] <- y 
    else {
      xgamma[i,j] <- xt
      kg[j] <- kg[j]+1 #y is rejected
    }
    vog[[j]]=xgamma[i,j]
    vog_x[[j]]=xgamma[i,j]
  }
  
}

###lambda
for (i in 1:21) {
  gamma[i]=mean(xgamma[1000:10000,i])
}
sum(gamma)

lamda_hat=mean(rgamma(100000,21+1,sum(gamma)+1))

##print
write.csv(c(kk,xbeta1,xbeta2),file.choose())
write.csv(kg,xgamma,file.choose())



