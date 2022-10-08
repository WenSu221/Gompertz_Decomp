
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
