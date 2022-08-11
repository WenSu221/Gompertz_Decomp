
###
###
###

library(DEoptim)

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
