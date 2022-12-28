
###
###
###

library(DEoptim)

# Arguments:

# Dx and Nx are matrix of the death counts and exposure
# Rows are ages and columns years
# The column should be named with the year, for example as "1950"
# age.start and age.end are the first and last age consider for the 
# decomposition (example: age.start=30 and age.end=100)



#### Gompertz-Makeham #####

Decomp.Makeham<- function(Dx,Nx, age.start, age.end){
  
  
  age<-c(age.start: age.end)
  
  
  ##names of the periods
  
  year1<- colnames(Dx[,1:(ncol(Dx)-1)])
  year2<- colnames(Dx[,2:ncol(Dx)])
  splitted <- (t(sapply(year2, function(x) substring(x, first=3, last=4))))
  period<- c(1:length(year1))
  for (i in 1:length(period)){
    period[i]<- paste(year1[i], splitted[i], sep="-")
  }
  
  
  ## Gompertz-Makeham, Poisson log-likelihood
  
  gomp.make <- function(c,M,b2,x){
    p1 <-  c + (b2*exp(b2*(x-M)))
    return(p1)
  }
  ll.poisson.make <- function(theta, Dx, Nx, x){
    c<- theta[1]
    M <- theta[2]
    b2 <- theta[3]
    out <-  -sum(Dx*log(gomp.make(x=x, c=c, M=M, b2=b2))-
                   gomp.make(x=x, c=c, M=M, b2=b2)*Nx)
    return(out)
    
  }
  par<- matrix(NA, nrow=ncol(Dx), ncol=3)
  for(i in 1:nrow(par)){
    par[i,] <- DEoptim(fn=ll.poisson.make, lower=c( 0, 50,  0.04),  
                       upper=c(0.1, 110,  4), 
                       Dx=Dx[,i], 
                       Nx= Nx[,i],x=age, 
                       control=DEoptim.control(trace=300))$optim$bestmem
  }
  rownames(par)<-colnames(Dx)
  colnames(par)<- c("c" ,"M", "B")
  
  
  ## Makeham term (c) contribution function (variability effect)
  
  c.ex.fun<-function(c1, c2, B1, B2, M1, M2){
    f2<-function(x){
      # Survival function
      lx1<- exp(-(c1*x) - (exp(-B1*M1)*(exp(B1*x)-1)))
      lx2<- exp(-(c2*x) - (exp(-B2*M2)*(exp(B2*x)-1)))
      # c dot
      c.dot<- log((c2/c1))*(c1*c2)^(1/2)
      #contributions
      cont<- -c.dot*x*(lx1*((lx2/lx1)^(1/2)))
      return(cont)
    }
    return(f2)
  }
  int.c <- function (c1, c2, B1, B2, M1, M2) { 
    integrate(c.ex.fun(c1, c2, B1, B2, M1, M2),  
              lower=age.start, upper=age.end)$value 
  }
  c.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(c.ex)){
    c.ex[i]<- int.c(c1=par[i,1], c2=par[(i+1),1],
                    B1=par[i,3], B2=par[(i+1),3], 
                    M1=par[i,2], M2=par[(i+1),2])
  }
  names(c.ex)<-period
  
  
  ## Beta contribution function (variability effect)
  
  B.ex.fun<-function(c1, c2, B1, B2, M1, M2){
    f2<-function(x){
      # Survival function
      lx1<- exp(-(c1*x) - (exp(-B1*M1)*(exp(B1*x)-1)))
      lx2<- exp(-(c2*x) - (exp(-B2*M2)*(exp(B2*x)-1)))
      # Cumulative hazard function for the Gompertz
      Hx1<-exp(-B1*M1)*(exp(B1*x)-1)
      Hx2<-exp(-B2*M2)*(exp(B2*x)-1)
      # Beta dot
      B.dot<- log(B2/B1)*(B1*B2)^(1/2)
      #contributions
      cont<- -B.dot*
        (lx1*((lx2/lx1)^(1/2)))*
        (((Hx1*(x-M1))+(x*exp(-B1*M1)))*
           abs(((Hx2*(x-M2))+(x*exp(-B2*M2)))/
                 ((Hx1*(x-M1))+(x*exp(-B1*M1))))^(1/2))
      return(cont)
    }
    return(f2)
  }
  int.B <- function ( c1, c2, B1, B2, M1, M2) { 
    integrate(B.ex.fun(c1, c2, B1, B2, M1, M2),  
              lower=age.start, upper=age.end)$value 
  }
  B.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(B.ex)){
    B.ex[i]<- int.B(c1=par[i,1], c2=par[(i+1),1],
                    B1=par[i,3], B2=par[(i+1),3], 
                    M1=par[i,2], M2=par[(i+1),2])
  }
  names(B.ex)<-period
  
  
  ##M contribution function (shifting effect)
  
  M.ex.fun<-function(c1, c2, B1, B2, M1, M2){
    f2<-function(x){
      # Survival
      lx1<- exp(-(c1*x) - (exp(-B1*M1)*(exp(B1*x)-1)))
      lx2<- exp(-(c2*x) - (exp(-B2*M2)*(exp(B2*x)-1)))
      # Cumulative hazard, Gompertz
      Hx1<-exp(-B1*M1)*(exp(B1*x)-1)
      Hx2<-exp(-B2*M2)*(exp(B2*x)-1)
      # M. dot
      M.dot<- log(M2/M1)*(M1*M2)^(1/2)
      #contributions
      cont<- M.dot*
        (lx1*((lx2/lx1)^(1/2)))*
        (Hx1*(Hx2/Hx1)^(1/2))*
        (B1*(B2/B1)^(1/2))
      return(cont)
    }
    return(f2)
  }
  int.M <- function (c1, c2, B1, B2, M1, M2) { 
    integrate(M.ex.fun(c1, c2, B1, B2, M1, M2),  
              lower=age.start, upper=age.end)$value 
  }
  M.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(M.ex)){
    M.ex[i]<- int.M(c1=par[i,1], c2=par[(i+1),1],
                    B1=par[i,3], B2=par[(i+1),3], 
                    M1=par[i,2], M2=par[(i+1),2])
  }
  names(M.ex)<-period
  
  
  ##Life expectancy
  
  lx.fun<- function(c,B,M){
    f2<- function(x){
      exp(-(c*x) - (exp(-B*M)*(exp(B*x)-1)))
    }
    return(f2)
  }
  int2 <- function (c,B,M) { 
    integrate(lx.fun(c,B,M),  lower=age.start, 
              upper=age.end)$value 
  }
  ex<- c(1:nrow(par))
  for(i in 1:length(ex)){
    ex[i]<- int2( c=par[i,1],  B=par[i,3], M=par[i,2])
  }
  names(ex)<-colnames(Dx)
  
  
  ##Changes in life expectancy
  
  diff.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(diff.ex)){
    diff.ex[i]<- ex[1+i]-ex[i]
  }
  names(diff.ex)<-period
  
  
  ##Values returned
  
  results<- list(par=par, deltac=c.ex, 
                 deltaB=B.ex, deltaM=M.ex, 
                 deltaex=diff.ex, ex=ex)
  return<- results
}



#### Siler #####

Decomp.Siler<- function(Dx,Nx, age.start, age.end){
  
  age=c(age.start:age.end)
  
  
  ##names of the periods
  
  year1<- colnames(Dx[,1:(ncol(Dx)-1)])
  year2<- colnames(Dx[,2:ncol(Dx)])
  splitted <- (t(sapply(year2, function(x) substring(x, first=3, last=4))))
  period<- c(1:length(year1))
  for (i in 1:length(period)){
    period[i]<- paste(year1[i], splitted[i], sep="-")
  }
  
  
  ## Siler, Poisson log-likelihood
  
  siler <- function(a1,b1,c,M,b2,x){
    p1 <- (a1*exp(-b1*x))+ c + (b2*exp(b2*(x-M)))
    return(p1)
  }
  ll.poisson.siler <- function(theta, Dx, Nx, x){
    a1<-theta[1]
    b1<- theta[2]
    c<- theta[3]
    M <- theta[4]
    b2 <- theta[5]
    out <-  -sum(Dx*log(siler(x=x, a1=a1, b1=b1, c=c, M=M, b2=b2))-
                   siler(x=x, a1=a1, b1=b1, c=c, M=M, b2=b2)*Nx,na.rm = T)
    return(out)
    
  }
  par<- matrix(NA, nrow=ncol(Dx), ncol=5)
  for(i in 1:nrow(par)){
    par[i,] <- DEoptim(fn=ll.poisson.siler, lower=c(0.000001, 0,  0, 50,  0.04),  
                       upper=c(0.7, 6,  0.1, 110,  1), 
                       Dx=Dx[,i], 
                       Nx= Nx[,i],x=age, 
                       control=DEoptim.control(trace=300))$optim$bestmem
  }
  rownames(par)<-colnames(Dx)
  colnames(par)<- c("a", "b" ,"c" ,"M", "B")
  
  
  ##a contribution function (variability effect)
  
  a.ex.fun<-function(a1,a2,b1,b2, c1, c2, B1, B2, M1, M2){
    f2<-function(x){
      #survivalfunction
      lx1<- exp((a1/b1*(exp(-b1*x)-1)) -(c1*x) - (exp(-B1*M1)*(exp(B1*x)-1)))
      lx2<- exp((a2/b2*(exp(-b2*x)-1)) -(c2*x) - (exp(-B2*M2)*(exp(B2*x)-1)))
      # a dot
      a.dot<- log(a2/a1)*(a1*a2)^(1/2)
      #contributions
      cont<- -a.dot *
        (lx1*((lx2/lx1)^(1/2)))*
        (-1/b1*((b1/b2)^(1/2)))*
        ((exp(-b1*x) - 1)*abs((exp(-b2*x) - 1)/(exp(-b1*x) - 1))^(1/2))
      return(cont)
    }
    return(f2)
  }
  int.a <- function (a1,a2,b1,b2, c1, c2, B1, B2, M1, M2) { 
    integrate(a.ex.fun(a1,a2,b1,b2, c1, c2, B1, B2, M1, M2),  lower=0, upper=110)$value 
  }
  a.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(a.ex)){
    a.ex[i]<- int.a(a1=par[i,1], a2=par[(i+1),1],
                    b1=par[i,2], b2=par[(i+1),2],
                    c1=par[i,3], c2=par[(i+1),3],
                    B1=par[i,5], B2=par[(i+1),5], 
                    M1=par[i,4], M2=par[(i+1),4])
  }
  names(a.ex)<-period
  
  
  ##b contribution function (variability effect)
  b.ex.fun<-function(a1,a2,b1,b2, c1, c2, B1, B2, M1, M2){
    f2<-function(x){
      #Survival function
      lx1<- exp((a1/b1*(exp(-b1*x)-1)) -(c1*x) - (exp(-B1*M1)*(exp(B1*x)-1)))
      lx2<- exp((a2/b2*(exp(-b2*x)-1)) -(c2*x) - (exp(-B2*M2)*(exp(B2*x)-1)))
      # b dot
      b.dot<- log(b2/b1)*(b1*b2)^(1/2)
      #contributions
      cont<- b.dot*
        (lx1*((lx2/lx1)^(1/2)))*
        (-a1*(-a2/-a1)^(1/2))*
        (1/b1*((1/b2)/(1/b1))^(1/2))*
        (((x*exp(-b1*x))+(exp(-b1*x)/b1)-(1/b1))*
           abs(((x*exp(-b2*x))+(exp(-b2*x)/b2)-(1/b2))/
                 ((x*exp(-b1*x))+(exp(-b1*x)/b1)-(1/b1)))^(1/2))
      return(cont)
    }
    return(f2)
  }
  int.b <- function (a1,a2,b1,b2, c1, c2, B1, B2, M1, M2) { 
    integrate(b.ex.fun(a1,a2,b1,b2, c1, c2, B1, B2, M1, M2),  lower=0, upper=110)$value 
  }
  b.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(b.ex)){
    b.ex[i]<- int.b(a1=par[i,1], a2=par[(i+1),1],
                    b1=par[i,2], b2=par[(i+1),2],
                    c1=par[i,3], c2=par[(i+1),3],
                    B1=par[i,5], B2=par[(i+1),5], 
                    M1=par[i,4], M2=par[(i+1),4])
  }
  names(b.ex)<-period
  
  
  ## Makeham term (c) contribution function (variability effect)
  
  c.ex.fun<-function(a1,a2,b1,b2, c1, c2, B1, B2, M1, M2){
    f2<-function(x){
      # Survival function
      lx1<- exp((a1/b1*(exp(-b1*x)-1)) -(c1*x) - (exp(-B1*M1)*(exp(B1*x)-1)))
      lx2<- exp((a2/b2*(exp(-b2*x)-1)) -(c2*x) - (exp(-B2*M2)*(exp(B2*x)-1)))
      # c dot
      c.dot<- log(c2/c1)*(c1*c2)^(1/2)
      #contributions
      cont<- -c.dot*x*(lx1*((lx2/lx1)^(1/2)))
      return(cont)
    }
    return(f2)
  }
  int.c <- function (a1,a2,b1,b2, c1, c2, B1, B2, M1, M2) { 
    integrate(c.ex.fun(a1,a2,b1,b2, c1, c2, B1, B2, M1, M2),  lower=0, upper=110)$value 
  }
  c.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(c.ex)){
    c.ex[i]<- int.c(a1=par[i,1], a2=par[(i+1),1],
                    b1=par[i,2], b2=par[(i+1),2],
                    c1=par[i,3], c2=par[(i+1),3],
                    B1=par[i,5], B2=par[(i+1),5], 
                    M1=par[i,4], M2=par[(i+1),4])
  }
  names(c.ex)<-period
  
  
  ## Beta contribution function (variability effect)
  
  B.ex.fun<-function(a1,a2,b1,b2, c1, c2, B1, B2, M1, M2){
    f2<-function(x){
      # Survival function
      lx1<- exp((a1/b1*(exp(-b1*x)-1)) -(c1*x) - (exp(-B1*M1)*(exp(B1*x)-1)))
      lx2<- exp((a2/b2*(exp(-b2*x)-1)) -(c2*x) - (exp(-B2*M2)*(exp(B2*x)-1)))
      # Cumulative hazard function for the Gompertz
      Hx1<-exp(-B1*M1)*(exp(B1*x)-1)
      Hx2<-exp(-B2*M2)*(exp(B2*x)-1)
      # Beta dot
      B.dot<- log(B2/B1)*(B1*B2)^(1/2)
      #contributions
      cont<- -B.dot*
        (lx1*((lx2/lx1)^(1/2)))*
        (((Hx1*(x-M1))+(x*exp(-B1*M1)))*
           abs(((Hx2*(x-M2))+(x*exp(-B2*M2)))/
                 ((Hx1*(x-M1))+(x*exp(-B1*M1))))^(1/2))
      return(cont)
    }
    return(f2)
  }
  int.B <- function (a1,a2,b1,b2, c1, c2, B1, B2, M1, M2) { 
    integrate(B.ex.fun(a1,a2,b1,b2, c1, c2, B1, B2, M1, M2),  lower=0, upper=110)$value 
  }
  B.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(B.ex)){
    B.ex[i]<- int.B(a1=par[i,1], a2=par[(i+1),1],
                    b1=par[i,2], b2=par[(i+1),2],
                    c1=par[i,3], c2=par[(i+1),3],
                    B1=par[i,5], B2=par[(i+1),5], 
                    M1=par[i,4], M2=par[(i+1),4])
  }
  names(B.ex)<-period
  
  
  ##Modal contribution function (shifting effect)
  
  M.ex.fun<-function(a1,a2,b1,b2, c1, c2, B1, B2, M1, M2){
    f2<-function(x){
      # Survival function
      lx1<- exp((a1/b1*(exp(-b1*x)-1)) -(c1*x) - (exp(-B1*M1)*(exp(B1*x)-1)))
      lx2<- exp((a2/b2*(exp(-b2*x)-1)) -(c2*x) - (exp(-B2*M2)*(exp(B2*x)-1)))
      # Cumulative hazard function for the Gompertz
      Hx1<-exp(-B1*M1)*(exp(B1*x)-1)
      Hx2<-exp(-B2*M2)*(exp(B2*x)-1)
      # M dot
      M.dot<- log(M2/M1)*(M1*M2)^(1/2)
      #contributions
      cont<- M.dot*
        (lx1*((lx2/lx1)^(1/2)))*
        (Hx1*(Hx2/Hx1)^(1/2))*
        (B1*(B2/B1)^(1/2))
      return(cont)
    }
    return(f2)
  }
  int.M <- function (a1,a2,b1,b2, c1, c2, B1, B2, M1, M2) { 
    integrate(M.ex.fun(a1,a2,b1,b2, c1, c2, B1, B2, M1, M2),  lower=0, upper=110)$value 
  }
  M.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(M.ex)){
    M.ex[i]<- int.M(a1=par[i,1], a2=par[(i+1),1],
                    b1=par[i,2], b2=par[(i+1),2],
                    c1=par[i,3], c2=par[(i+1),3],
                    B1=par[i,5], B2=par[(i+1),5], 
                    M1=par[i,4], M2=par[(i+1),4])
  }
  names(M.ex)<-period
  
  
  ##Life expectancy
  
  lx.fun<- function(a,b,c,B,M){
    f2<- function(x){
      exp((a/b*(exp(-b*x)-1)) -(c*x) - (exp(-B*M)*(exp(B*x)-1)))
    }
    return(f2)
  }
  int2 <- function (a,b,c,B,M) { 
    integrate(lx.fun(a,b,c,B,M),  lower=0, upper=110)$value 
  }
  ex<- c(1:nrow(par))
  for(i in 1:length(ex)){
    ex[i]<- int2( a=par[i,1], b=par[i,2], c=par[i,3],  B=par[i,5], M=par[i,4])
  }
  names(ex)<-colnames(Dx)
  
  
  ##Changes in life expectancy
  
  diff.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(diff.ex)){
    diff.ex[i]<- ex[1+i]-ex[i]
  }
  names(diff.ex)<-period
  
  
  ## Values returned
  
  results<- list(par=par, deltaa=a.ex, deltab=b.ex, deltac=c.ex, deltaB=B.ex, deltaM=M.ex, 
                 deltaex=diff.ex, ex=ex)
  return<- results
}
