#defenition of I&M-splines- Interval
Mspline=function(l,t,k,m=20){
  a=0
  R=LR[,2]
  b=max(R[is.finite(R)])
  kesi=numeric(m)
  for (i in 1:m) {
    kesi[i]=i*((b-a)/(m+1))
  }
  s <- numeric(m+2*k)
  for (i in 1:(k)) {
    s[i]=a
  }
  for (i in (k+1):(m+k)) {
    s[i] = kesi[i-k]
  }
  for (i in (m+k+1):(m+2*k)) {
    s[i]=b
  }
  print(s)
  if (k==1){
    if (s[l]<=t & t<s[l+1]){
      print(c(s[l+1],s[l]))
      sp = (1)/(s[l+1]-s[l])
    }
    else {sp = 0
    }
  }
  else {
    if (s[l]<=t & t<s[l+k]){
      print(c(s[l+k],s[l]))
      sp1 = ((t-s[l])*Mspline(l,t,k-1)+(s[l+k]-t)*Mspline(l+1,t,k-1))/
        (s[l+k]-s[l])
      print(k)
      sp = sp1*k/(k-1)
    }
    else {sp = 0}
  }
  sp
}

Ispline=function(l,t,k,m=20){
  a=0
  R=LR[,2]
  d=k-1
  sp=0
  b=max(R[is.finite(R)])
  kesi=numeric(m)
  for (i in 1:m) {
    kesi[i]=i*((b-a)/(m+1))
  }
  s <- numeric(m+2*k+2)
  for (i in 1:(k+1)) {
    s[i]=a
  }
  for (i in (k+2):(m+k+1)) {
    s[i] = kesi[i-k]
  }
  for (i in (m+k+2):(m+2*k+2)) {
    s[i]=b
  }
  for (i in 1:(m+2*k+1)) {
    if(s[i]<= t & t<s[i+1]) j=i
  }
  print(s)
  if(l>j) {sp=0}
  else if (l<j-k+1) {sp=1}
  else if (j-k+1 <= l & l<= j){
    for (m in l:j) {
      sp=sp+(s[m+k+1]-s[m])*Mspline(m,t,k+1)/(k+1)
    }
  }
  sp
}
