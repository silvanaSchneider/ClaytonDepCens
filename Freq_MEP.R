
###---------------------------------------------------------------------------------------------###
# Functions for the Clayton copula model with piecewise exponential marginal distribution,        #              # 
# with inference via maximum likelihood.                                                          #
###---------------------------------------------------------------------------------------------###

rm(list=ls(all=TRUE)) 
wd <- getwd()

###---------------------------------------------------------------------------------------------###
# Functions to generate datasets from Clayton copula with piecewise exponential distribution      #
###---------------------------------------------------------------------------------------------###

library(gamlss);library(gamlss.dist);library(PWEALL)

###---------------------------------------------------------------------------------------------###
# Generate datasets
###---------------------------------------------------------------------------------------------###

k=500  # number of data sets
n=200  # sample size

### Parameter Values
alfaa=3

a_t=2
a_c=2

l_t=1.5   
l_c=1.5

CA=1.2

beta_t=c(-1,1.3)
beta_c=c(-1,1.3)

datasets=list()


for(j in 1:k){
  datasets[[j]]=list()
  
  set.seed(j*7+20)
  
  x_t=matrix(c(rbinom(n,1,0.5),rnorm(n,0,1)),ncol = 2,byrow = F)  
  x_c=matrix(c(rbinom(n,1,0.5),rnorm(n,0,1)),ncol = 2,byrow = F)
  
  eta_t=exp(x_t%*%beta_t)
  eta_c=exp(x_c%*%beta_c)
  
  U <- cbind(runif(n),runif(n)) # Uniform variables used in Clayton copula
  U[,2] <- (((U[,1]^(-alfaa))*((U[,2]^(-alfaa/(1+alfaa)))-1))+1)^(-1/alfaa)
  
  ## Generating Clayton copula Weibull model data
  t=vector(length = n) # lifetimes
  for (i in 1:n) { t[i]=(log(1-U[i,1])/(-l_t*eta_t[i]))^(1/a_t) }
  
  c=vector(length = n) # censoring times
  for (i in 1:n) { c[i]=(log(1-U[i,2])/(-l_c*eta_c[i]))^(1/a_c)}
  
  data=mapply(min,t[1:n],c[1:n],CA) ## take minimum between times
  ind=(t[1:n]<=c[1:n] & t[1:n]<=CA )*1 ## event indicator: 1 if is failure; 0 if is censored
  rho=(t[1:n]<=CA |c[1:n] <= CA )*1
  
  datasets[[j]][[1]]=cbind(t,c)
  datasets[[j]][[2]]=cbind(data,ind,rho)
  datasets[[j]][[3]]=x_t
  datasets[[j]][[4]]=x_c
  
}

###---------------------------------------------------------------------------------------------###
# Model
###---------------------------------------------------------------------------------------------###

#---------------------------------------------------------------------------------------------
# Functions used in the piecewise Exponential distribution

time.grid <- function(time, event, n.int=NULL)
{
  o <- order(time)  
  time <- time[o]    
  event <- event[o]
  time.aux <- unique(time[event==1])
  if(is.null(n.int))
  {
    n.int <- length(time.aux)
  }
  
  m <- length(time.aux)
  if(n.int > m)
  {
    a <- c(0,unique(time[event==1]))
    a[length(a)] <- Inf
  }
  else
  {
    b <- min(m,n.int)
    k1 <- trunc(m/b)
    r <- m-b*k1
    k2 <- k1+1
    idf1 <- seq(k1,(b-r)*k1, k1)
    idf2 <- sort(seq(m,max(idf1),-k2))
    idf <- unique(c(idf1,idf2))
    a_inf <- c(0,time.aux[idf])
    a_inf[length(a_inf)] <- Inf
    a_s_inf  <- c(0,time.aux[idf])
  }
  saida <- list(a_inf,a_s_inf)
  
  return(saida)
}


#---------------------------------------------------------------------------------------------
# Hazard functions 
h_t=function(lambda.t,beta_te) lambda.t[id.falha]*exp(x_t%*%beta_te)
h_c=function(lambda.c,beta_ce) lambda.c[id.cens]*exp(x_c%*%beta_ce)

H_t=function(lambda.t,beta_te) exp(x_t%*%beta_te)*(xi.falha%*%lambda.t)
H_c=function(lambda.c,beta_ce) exp(x_c%*%beta_ce)*(xi.cens%*%lambda.c) 


# Clayton copula              
C_a= function(H_t,H_c) 1-(1-exp(-H_t))-(1-exp(-H_c))+(( (1-exp(-H_t))^(-alfa)+(1-exp(-H_c))^(-alfa)-1)^(-1/alfa)) ##equacao da copula #1-u-v+C(u,v), onde u e v sao Funcoes de ditribuicao acumuladas, essa formula??o compoe uma copula de sobrevivencia

dC_t=function(H_t,H_c)((1-exp(-H_t))^(-alfa-1))*( (1-exp(-H_t))^(-alfa)+(1-exp(-H_c))^(-alfa)-1)^(-1/alfa-1) -1 ##derivada da copula em rela??o a S_t #derivada (1-u-v+C(u,v)) em rela??o a u

dC_c=function(H_t,H_c)((1-exp(-H_c))^(-alfa-1))*( (1-exp(-H_t))^(-alfa)+(1-exp(-H_c))^(-alfa)-1)^(-1/alfa-1) -1   ##derivada da copula em rela??o a S_c #derivada (1-u-v+C(u,v)) em rela??o a v

#---------------------------------------------------------------------------------------------
# Likelihood function 
veros=function(theta){
  
  lambda.t=exp(theta[1:bmax_t])
  
  lambda.c=exp(theta[(bmax_t+1):(bmax_t+bmax_c)])
  
  beta_te=theta[(bmax_t+bmax_c+1):(bmax_t+bmax_c+p)]
  beta_ce=theta[(bmax_t+bmax_c+p+1):(bmax_t+bmax_c+p+q)]
  
  logl= sum(   rho*(ind * ( log(h_t(lambda.t,beta_te))-H_t(lambda.t,beta_te)+log(-dC_t(H_t(lambda.t,beta_te),H_c(lambda.c,beta_ce))+0.0000001)  )
                    +(1-ind)* ( log(h_c(lambda.c,beta_ce))-H_c(lambda.c,beta_ce)+log(-dC_c(H_t(lambda.t,beta_te),H_c(lambda.c,beta_ce))+0.0000001)))+
                 (1-rho)* log(C_a(H_t(lambda.t,beta_te),H_c(lambda.c,beta_ce))+0.0000001))
  
  return(-logl)
}

###---------------------------------------------------------------------------------------------###
# Fit with alpha=3 (true value)
###---------------------------------------------------------------------------------------------###

time=datasets[[1]][[2]][,1] 
ind=datasets[[1]][[2]][,2]
rho=datasets[[1]][[2]][,3]
x_t=datasets[[1]][[3]]
x_c=datasets[[1]][[4]]
p = ncol(x_t)
q = ncol(x_c)
n = length(time)

param_real=c(log(c(a_t,l_t,a_c,l_c)),beta_t,beta_c)   ##vector of parameter values, generated data Weibull

bmax_t  <- 5
bmax_c  <- 5

alfa <- 3
theta <- rep(0.1,p+q+bmax_t+bmax_c) # inicial values

est=matrix(NA,nrow = k,ncol = (p+q+bmax_t+bmax_c))
hessiana=matrix(NA,nrow = k,ncol = (p+q+bmax_t+bmax_c))
loglik=matrix(NA,nrow = k,ncol = 1)

for(l in 1:k){
  
  time=datasets[[l]][[2]][,1] 
  ind=datasets[[l]][[2]][,2]
  rho=datasets[[l]][[2]][,3]
  x_t=datasets[[l]][[3]]
  x_c=datasets[[l]][[4]]
  p = ncol(x_t)
  q = ncol(x_c)
  n = length(time)
  
  delta.t <- ind
  delta.c <- (1-ind)*rho
  
  a <- time.grid(time,delta.t,bmax_t)
  a <- a[[1]]
  s <- time.grid(time,delta.c,bmax_c) 
  s <- s[[1]]
  b <- length(a)-1
  c <- length(s)-1
  a[length(a)] <- max(time)
  s[length(s)] <- max(time)
  
  xi.falha <- matrix(0,n,b)
  xi.cens <- matrix(0,n,c)
  t <- time
  
  for(i in 1:n){                         
    for(j in 1:b){
      xi.falha[i,j] <- (min(t[i], a[j+1])-a[j])*((t[i]-a[j])>0)
    }
  }    
  for(i in 1:n){                        
    for(j in 1:c){
      xi.cens[i,j] <- (min(t[i], s[j+1])-s[j])*((t[i]-s[j])>0)
    }
  }                                                 
  
  id.falha <- as.numeric(cut(t,a))
  id.cens  <- as.numeric(cut(t,s))
  
  tt= optim(rep(0.01, (p+q+bmax_t+bmax_c)),veros)
  tt= optim(tt$par,veros)
  tt= optim(tt$par,veros,hessian = T)
  hessiana[l,]=sqrt(diag(solve(tt$hessian)))
  est[l,]=c(exp(tt$par[1:(bmax_t+bmax_c)]),tt$par[(bmax_t+bmax_c+1):(bmax_t+bmax_c+p+q)])
  loglik[l,] = tt$value
  print(l)
}


# saving the estimated results 
write.table(est,file='Est_MEP_alpha3.txt', col.names = F, row.names = F)   # estimates saved in sequence c(lambda.t, lambda.c ,beta_t,beta_c)
write.table(hessiana,file='EP_MEP_alpha3.txt', col.names = F, row.names = F)  # standard error saved in sequence c(lambda.t, lambda.c ,beta_t,beta_c)
write.table(loglik,file='Loglik_EP_alpha3.txt', col.names = F, row.names = F)

###---------------------------------------------------------------------------------------------###
# Analysis of outputs, assuming alha=3

est1      <- est[,(bmax_t+bmax_c+1):(bmax_t+bmax_c+p+q)]
ErroP     <- hessiana[,(bmax_t+bmax_c+1):(bmax_t+bmax_c+p+q)]

para_real <- c(param_real[5:6],param_real[7:8]) # only evaluate beta_t and beta_c as data is generated from Weibull

RV <- (est1-matrix(rep(para_real,k),nrow = k,byrow =T))*100/matrix(rep(para_real,k),nrow = k,byrow =T)  
boxplot(RV)

ampitude <-  (est1+1.96*ErroP)-(est1-1.96*ErroP)   

#coverage probability
R <- nrow(est1)
lim_inf_r <- matrix(NA,nrow(est1),4)
lim_sup_r <- matrix(NA,nrow(est1),4)

for (i in 1:4){
  lim_inf_r[,i] <- est1[,i] -1.96*sd(est1[,i])
  lim_sup_r[,i] <- est1[,i] +1.96*sd(est1[,i])
}
Prob_cobertura <- NULL
for (i in 1:4){
  Prob_cobertura[i] <- (sum(lim_sup_r[,i] >= para_real[i] & lim_inf_r[,i] <= para_real[i]))/R
}


# Table
stat<-matrix(rep(NA,7*4),ncol=7,nrow=4)
stat<-as.data.frame(stat)
row.names(stat) <- c("$\\beta^{T}_{1}$","$\\beta^{T}_{2}$","$\\beta^{C}_{1}$","$\\beta^{C}_{2}$")  #$\\alpha$: fixo
colnames(stat)<- c( "Real", "Estimate", "SE", "SD", "RV%", "Amplitude", "PC")                      #SE:standard error; RV\% vicio relativo*100; Amplitude:amplitude do IC
stat[,1] <- para_real
stat[,2] <- paste(round(apply(est1[,],2,mean),3))
stat[,3] <- paste(round(apply(ErroP[,],2,mean),3)) # average of 500 standard errors
stat[,4] <- paste(round(apply(est1[,],2,sd),3))  # standard deviation
stat[,5] <- paste(round(apply(RV[,],2,mean),3)) # RV\% relavite bias
stat[,6] <- paste(round(apply(ampitude[,],2,mean),3)) # Amplitude
stat[,7] <- round(Prob_cobertura,3) # coverage probability

write.table(stat,file='Tab_MEP_alpha3.txt') 
