
####
#### Adding confidence interval for Siler model
####

library(data.table)
library(DEoptim)

#### Functions ####

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

#### Data fit ####

Dx <- fread("data/USA.Deaths_1x1.txt")[,list(Year,Age,Female)]

Dx <- matrix(data = Dx$Female,nrow = 111,ncol = length(unique(Dx$Year)),
             dimnames = list(unique(Dx$Age),unique(Dx$Year)))

Ex <- fread("data/USA.Exposures_1x1.txt")[,list(Year,Age,Female)]

Ex <- matrix(data = Ex$Female,nrow = 111,ncol = length(unique(Ex$Year)),
             dimnames = list(unique(Ex$Age),unique(Ex$Year)))

par_mat<- matrix(NA, nrow=ncol(Dx), ncol=5)

age <- 0:110 

for(i in 1:nrow(par_mat)){
  par_mat[i,] <- DEoptim(fn=ll.poisson.siler, lower=c(0.000001, 0,  0, 50,  0.04),  
                     upper=c(0.7, 6,  0.1, 110,  1), 
                     Dx=Dx[,i], 
                     Nx= Ex[,i],
                     x=age, 
                     control=DEoptim.control(NP=100,
                                             trace=300))$optim$bestmem
}

#### Bootstrap the mortality rate ####

age <- 0:110

m <- nrow(Dx)
Ns <- 100

qx <- 1-exp((Dx[,1]/Ex[,1])*-1)

boot_mat <- suppressWarnings(matrix(rbinom(m*Ns,
                                    as.integer(Ex[,1]),
                                    qx),
                             nrow=m,ncol=Ns))

table_dth <- data.table(parameter = 1:5)

for(i in 1:ncol(boot_mat)){
  vec <- DEoptim(fn=ll.poisson.siler, 
                         lower=c(0.000001, 0,  0, 50,  0.04),  
                         upper=c(0.7, 6,  0.1, 110,  1), 
                         Dx=boot_mat[,i], 
                         Nx= Ex[,1],
                         x=age, 
                         control=DEoptim.control(NP=100,
                                                 trace=300))$optim$bestmem
  table_dth[,paste0(i):=vec]
}


hist(unlist(table_dth[4,-1]))
