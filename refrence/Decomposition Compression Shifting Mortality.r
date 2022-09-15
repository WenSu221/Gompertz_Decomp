###		WRITTEN BY MARIE-PIER BERGERON-BOUCHER

##This code are for assessing differences between 2 consecutive years only

##### the following libraries will be needed!!!!

library(DEoptim)
library(RColorBrewer)
library(lattice)
library(dichromat)
library(fds)


##################   DATA     #####################################


##Extracting data from the Human Mortality Database (HMD), Example for Sweden
#Enter username and password for the HMD in the function hmdcountry()


sweden<- hmdcountry("SWE", sex="Female", 
                    username="u6897805@anu.edu.au", 
                    password="jYHy!m!6i5Ae6!8")
summary(sweden)

Dx<-sweden$Deathcount$y
Ex<-sweden$Exposure$y



#################################################################

## Gompertz, Gompertz-Makeham and Siler decomposition of life expectancy functions

# Description: The functions Decomp.Gompertz, Decomp.Makeham and Decomp.Siler decompose life expectancy
# at birth due to changes in the models parameters.

# Arguments:

# Dx and Nx are matrix of the death counts and exposure
# Rows are ages and columns years
# The column should be named with the year, for example as "1950"
# age.start and age.end are the first and last age consider for the decomposition (example: age.start=30 and age.end=100)

#################################################################



####################	Functions    ###########################3

##------------------------------------------------------------------------------------------##

## Gompertz decomposition of life expectancy


Decomp.Gompertz<- function(Dx,Nx, age.start, age.end){
  

  age<-c(age.start:age.end)
  
  
  ##names of the periods

  year1<- colnames(Dx[,1:(ncol(Dx)-1)])
  year2<- colnames(Dx[,2:ncol(Dx)])
  splitted <- (t(sapply(year2, function(x) substring(x, first=3, last=4))))
  period<- c(1:length(year1))
  for (i in 1:length(period)){
    period[i]<- paste(year1[i], splitted[i], sep="-")
  }
  
  
  ## Gompertz, Poisson log-likelihood
  
  gompertz <- function(M,b,x){
    p1 <- b*exp(b*(x-M))
    return(p1)
  }
  
  ll.poisson.gompertz <- function(theta, Dx, Nx, x){
    M <- theta[1]
    b <- theta[2]
    out <-  -sum(Dx*log(gompertz(x=x, M=M, b=b))-gompertz(x=x, M=M,b=b)*Nx)
    return(out)
  }
  
  par <- matrix(NA, nrow=ncol(Dx), ncol=2)
  
  for(i in 1:nrow(par)){
    par[i,] <- DEoptim(fn=ll.poisson.gompertz, lower=c(50,  0.01),  
                              upper=c(110,  4), Dx=Dx[,i], 
                              Nx=Nx[,i],x=age, 
                              control=DEoptim.control(trace=300))$optim$bestmem
  }
  
  rownames(par)<-colnames(Dx)
  colnames(par)<- c("M", "B")
  
  
  ## Beta contribution function (variability effect)
  
  B.ex.fun<- function(b1, b2, M1, M2){
    f2<-function(x){
      #Survival function
      lx1<-exp(-exp(-b1*M1)*(exp(b1*x)-1))
      lx2<-exp(-exp(-b2*M2)*(exp(b2*x)-1))
      #Cumulative hazard function
      Hx1<-exp(-b1*M1)*(exp(b1*x)-1)
      Hx2<-exp(-b2*M2)*(exp(b2*x)-1)
      #Beta dot
      B.dot<- log(b2/b1)*(b1*b2)^(1/2)
      #Beta contribution
      cont<- -B.dot*
        (lx1*((lx2/lx1)^(1/2)))*
        (((Hx1*(x-M1))+(x*exp(-b1*M1)))*
           abs(((Hx2*(x-M2))+(x*exp(-b2*M2)))/
                 ((Hx1*(x-M1))+(x*exp(-b1*M1))))^(1/2))
      return(cont)
    }
    return(f2)
  }
  int.B <- function (b1, b2, M1, M2) { 
    integrate(B.ex.fun(b1, b2, M1, M2),  lower=age.start, upper=age.end)$value 
  }
  B.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(B.ex)){
    B.ex[i]<- int.B(b1=par[i,2], b2=par[(i+1),2], 
                           M1=par[i,1], M2=par[(i+1),1])
  }
  names(B.ex)<-period
  
  
  ## M contribution function (shifting effect)
  
  M.ex.fun<- function(b1, b2, M1, M2){
    f2<-function(x){
      #Survival function
      lx1<-exp(-exp(-b1*M1)*(exp(b1*x)-1))
      lx2<-exp(-exp(-b2*M2)*(exp(b2*x)-1))
      #Cumulative hazard function
      Hx1<-exp(-b1*M1)*(exp(b1*x)-1)
      Hx2<-exp(-b2*M2)*(exp(b2*x)-1)
      #M dot
      M.dot<- log(M2/M1)*(M1*M2)^(1/2)
      #M contribution
      cont<- M.dot*
        (((lx2*lx1)^(1/2))*
           ((Hx2*Hx1)^(1/2))*
           ((b2*b1)^(1/2)))
      return(cont)
    }
    return(f2)
  }
  int.M <- function (b1, b2, M1, M2) { 
    integrate(M.ex.fun(b1, b2, M1, M2),  lower=age.start, upper=age.end)$value 
  }
  M.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(M.ex)){
    M.ex[i]<- int.M(b1=par[i,2], b2=par[(i+1),2], 
                           M1=par[i,1], M2=par[(i+1),1])
  }
  names(M.ex)<-period
  
 
  ##Life expectancy
 
  lx.fun<- function(b,M){
    f2<- function(x){
      exp(-exp(-b*M)*(exp(b*x)-1))
    }
    return(f2)
  }
  int2 <- function (b, M) { 
    integrate(lx.fun(b, M),  lower=age.start, upper=age.end)$value 
  }
  ex<- c(1:nrow(par))
  for(i in 1:length(ex)){
    ex[i]<- int2(par[i,2], par[i,1])
  }
  names(ex)<-colnames(Dx)
  
  
  ##Changes in life expectancy
  
  diff.ex<- c(1:(nrow(par)-1))
  for(i in 1:length(diff.ex)){
    diff.ex[i]<- ex[1+i]-ex[i]
  }
  names(diff.ex)<-period
 
  
  ##Returned values
  
  results<- list(par=par, deltaB=B.ex, deltaM=M.ex, deltaex=diff.ex, ex=ex)
  return<- results
}




##------------------------------------------------------------------------------------------##

## Gompertz-Makeham decomposition of life expectancy


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
    integrate(c.ex.fun(c1, c2, B1, B2, M1, M2),  lower=age.start, upper=age.end)$value 
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
    integrate(B.ex.fun(c1, c2, B1, B2, M1, M2),  lower=age.start, upper=age.end)$value 
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
    integrate(M.ex.fun(c1, c2, B1, B2, M1, M2),  lower=age.start, upper=age.end)$value 
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
    integrate(lx.fun(c,B,M),  lower=age.start, upper=age.end)$value 
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
  
  results<- list(par=par, deltac=c.ex, deltaB=B.ex, deltaM=M.ex, deltaex=diff.ex, ex=ex)
  return<- results
}




##------------------------------------------------------------------------------------------##

## Siler decomposition of life expectancy


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
                   siler(x=x, a1=a1, b1=b1, c=c, M=M, b2=b2)*Nx)
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




#------------------------------------------------------------------------------------------#



##################   RESULTS    #####################################

##Value returned: 

#deltaa: contributions to the change in life expectancy from the parameter a (Siler).
#deltab: contributions to the change in life expectancy from the parameter b (Siler).
#deltac: contributions to the change in life expectancy from the parameter c (Gompertz-Makeham and Siler).
#deltaB: contributions to the change in life expectancy from the parameter beta (Gompertz, Gompertz-Makeham and Siler).
#deltaM: contributions to the change in life expectancy from the parameter M (Gompertz, Gompertz-Makeham and Siler).
#delatex: increase in life expectancy between 2 consecutive years
#par: parameter of the selected model
#ex: life expectancy at birth

#################################################################


##Preparing the data for Gompertz and Gompertz-Makeham

age.interval <- c(30,110)
age <- (age.interval[1]:age.interval[2])
year.interval<-c(1900,2010)

death.counts <- Dx[(age.interval[1]+1):(age.interval[2]+1), as.character(c(year.interval[1]:year.interval[2]))]
expo.counts <- Ex[(age.interval[1]+1):(age.interval[2]+1), as.character(c(year.interval[1]:year.interval[2]))]


##Preparing the data for Siler

age.interval.siler <- c(0,110)
age.siler <- (age.interval.siler[1]:age.interval.siler[2])

death.counts.siler <- Dx[(age.interval.siler[1]+1):(age.interval.siler[2]+1), as.character(c(year.interval[1]:year.interval[2]))]
expo.counts.siler <- Ex[(age.interval.siler[1]+1):(age.interval.siler[2]+1), as.character(c(year.interval[1]:year.interval[2]))]


## Run codes

swe.gomp<- Decomp.Gompertz(Dx=death.counts, Nx=expo.counts, 
                        age.start=age.interval[1], age.end=age.interval[2])
swe.make<- Decomp.Makeham(Dx=death.counts, Nx=expo.counts, 
                           age.start=age.interval[1], age.end=age.interval[2])
swe.siler<- Decomp.Siler(Dx=death.counts.siler, Nx=expo.counts.siler, 
                          age.start=age.interval.siler[1], age.end=age.interval.siler[2])


## Results 1 year-period

swe.gomp$deltaB["2005-06"]+swe.gomp$deltaM["2005-06"]
swe.gomp$deltae["2005-06"]

swe.make$deltac["2005-06"]+swe.make$deltaB["2005-06"]+swe.make$deltaM["2005-06"]
swe.make$deltae["2005-06"]

swe.siler$deltaa["2005-06"]+swe.siler$deltab["2005-06"]+swe.siler$deltac["2005-06"]+swe.siler$deltaB["2005-06"]+swe.siler$deltaM["2005-06"]
swe.siler$deltae["2005-06"]


## Results 5 years period

v <- swe.gomp$deltaB
n <- 5
cutpoints <- seq( 1 , length( v ) , by = n )
categories <- findInterval( 1:length( v ) , cutpoints )

B5.ex.gomp<- tapply(swe.gomp$deltaB , categories , sum )
M5.ex.gomp<- tapply(swe.gomp$deltaM , categories , sum )
e5.ex.gomp<- tapply(swe.gomp$deltae , categories , sum )

c5.ex.make<- tapply(swe.make$deltac , categories , sum )
B5.ex.make<- tapply(swe.make$deltaB , categories , sum )
M5.ex.make<- tapply(swe.make$deltaM , categories , sum )
e5.ex.make<- tapply(swe.make$deltae , categories , sum )

a5.ex.siler<- tapply(swe.siler$deltaa , categories , sum )
b5.ex.siler<- tapply(swe.siler$deltab , categories , sum )
c5.ex.siler<- tapply(swe.siler$deltac , categories , sum )
B5.ex.siler<- tapply(swe.siler$deltaB , categories , sum )
M5.ex.siler<- tapply(swe.siler$deltaM , categories , sum )
e5.ex.siler<- tapply(swe.siler$deltae , categories , sum )



##Graph

year1<- seq(year.interval[1], year.interval[2]-5, 5)
year2<- c(seq(year.interval[1]+5, year.interval[2], 5))
splitted <- (t(sapply(year2, function(x) substring(x, first=3, last=4))))
period<- c(1:length(year1))
for (i in 1:length(period)){
  period[i]<- paste(year1[i], splitted[i], sep="-")
}

contributions.gomp<- cbind(B5.ex.gomp, M5.ex.gomp)
rownames(contributions.gomp)<-period

contributions.make<- cbind(c5.ex.make, B5.ex.make, M5.ex.make)
rownames(contributions.make)<-period

contributions.siler<- cbind(a5.ex.siler, b5.ex.siler, c5.ex.siler, B5.ex.siler, M5.ex.siler)
rownames(contributions.siler)<-period

mycol<-c(colorschemes$DarkRedtoBlue.18[c(1,3,5,7)],"#FD8D3C")


#Gompertz contributions to changes in life expectancy

barchart(as.matrix(contributions.gomp), stack=T, ylim=c(-1.5,3),
         col=mycol[c(1,5)], horizontal=F, 
         box.ratio=30, scales = list(x = list(rot = 45), cex=2,tick.number=10), xlab=list("Period", cex=2),  
         ylab=list("Contributions to change in life expectancy at age 30", cex=2),
         auto.key=list(columns=1, rectangles=T, x=0.7, y=0.2, 
                       text=c(expression(~ Delta ~ beta), 
                              expression(~ Delta ~ "M")), cex=3),
         par.settings = simpleTheme(col = mycol[c(1,5)]))


#Gompertz-Makeham contributions to changes in life expectancy

barchart(as.matrix(contributions.make), stack=T, ylim=c(-1.5,3),
         col=mycol[c(2,1,5)], horizontal=F, 
         box.ratio=30, scales = list(x = list(rot = 45), cex=2,tick.number=10), xlab=list("Period", cex=2),  
         ylab=list("Contributions to change in life expectancy at age 30", cex=2),
         auto.key=list(columns=1, rectangles=T, x=0.7, y=0.2, 
                       text=c(expression(~ Delta ~ "c"),
                              expression(~ Delta ~ beta), 
                              expression(~ Delta ~ "M")), cex=3),
         par.settings = simpleTheme(col = mycol[c(2,1,5)]))


#Siler contributions to changes in life expectancy

barchart(as.matrix(contributions.siler), stack=T, ylim=c(-1.5,5),
         col=mycol[c(4,3,2,1,5)], horizontal=F, 
         box.ratio=30, scales = list(x = list(rot = 45), cex=2,tick.number=10), xlab=list("Period", cex=2),  
         ylab=list("Contributions to change in life expectancy at age 0", cex=2),
         auto.key=list(columns=1, rectangles=T, x=0.8, y=0.9, 
                       text=c(expression(~ Delta ~ "a"),
                              expression(~ Delta ~ "b"),
                              expression(~ Delta ~ "c"),
                              expression(~ Delta ~ beta), 
                              expression(~ Delta ~ "M")), cex=3),
         par.settings = simpleTheme(col = mycol[c(4,3,2,1,5)]))


