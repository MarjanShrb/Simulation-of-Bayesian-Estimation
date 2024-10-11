x=seq(0.01,1,0.01)
y=numeric(100)
for(i in 1:100){
  y[i]= f(x[i])
  i=i+1
}
plot(x,y)
#algorithm
#f=function(x) if(x>0 & x<1) 6*x*(1-x) else return("NA")

L <- function(x, x1, y1, x2, y2) {
  tryCatch({
    ((y2 - y1) / (x2 - x1)) * (x - x1) + y1
  }, error = function(e) {
    return(NA)
  })
}
#define sn and lg: lg is log[f[x]]
h <- Vectorize(function(x) {
  i <- max(which(sn <= x))
  L1 <- L(x, sn[i - 1], lg[i - 1], sn[i], lg[i])
  L2 <- L(x, sn[i + 1], lg[i + 1], sn[i + 2], lg[i + 2])
  min(L1, L2, na.rm = TRUE)
})
hexp=Vectorize(function(x)exp(h(x)))
g=function(x)(1/integrate(hexp,g_low,g_up)$value)*hexp(x)


XA=c()
#i=2
sn=c(0,0.1,0.5)
kkk=0

while(length(XA)<1) {
  u=runif(1)
  X=uniroot(Vectorize(function(x)integrate(g,0,x)$value-u),c(0,10))$root
  U=runif(1)
  if (U<(f(X)/hexp(X))){
    XA=c(XA,X)
  } else{
    sn=c(sn,X)
    sn=sort(sn)
    kkk=kkk+1
  }
}
length(XA)




