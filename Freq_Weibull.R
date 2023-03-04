
###---------------------------------------------------------------------------------------------###
# Functions for the Clayton copula model with Weibull marginal distribution,                      # 
# with inference via maximum likelihood.                                                          #
###---------------------------------------------------------------------------------------------###

rm(list=ls(all=TRUE)) 
wd <- getwd()


###---------------------------------------------------------------------------------------------###
# Functions to generate datasets from Clayton copula with Weibull distribution for times          #
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
  t=vector(length = n) # lifetime
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
# Hazard functions 
h_t=function(a_t,l_t,beta_te) a_t*l_t*data^(a_t-1)*exp(covariaveis_t%*%beta_te)
h_c=function(a_c,l_c,beta_ce) a_c*l_c*data^(a_c-1)*exp(covariaveis_c%*%beta_ce)

H_t=function(a_t,l_t,beta_te) l_t*data^(a_t)*exp(covariaveis_t%*%beta_te)
H_c=function(a_c,l_c,beta_ce) l_c*data^(a_c)*exp(covariaveis_c%*%beta_ce)

# Clayton copula         
C_a= function(H_t,H_c) 1-(1-exp(-H_t))-(1-exp(-H_c))+(( (1-exp(-H_t))^(-alfa)+(1-exp(-H_c))^(-alfa)-1)^(-1/alfa)) ##equacao da copula #1-u-v+C(u,v), onde u e v sao Funcoes de ditribuicao acumuladas, essa formula??o compoe uma copula de sobrevivencia

dC_t=function(H_t,H_c)((1-exp(-H_t))^(-alfa-1))*( (1-exp(-H_t))^(-alfa)+(1-exp(-H_c))^(-alfa)-1)^(-1/alfa-1) -1 ##derivada da copula em rela??o a S_t #derivada (1-u-v+C(u,v)) em rela??o a u

dC_c=function(H_t,H_c)((1-exp(-H_c))^(-alfa-1))*( (1-exp(-H_t))^(-alfa)+(1-exp(-H_c))^(-alfa)-1)^(-1/alfa-1) -1   ##derivada da copula em rela??o a S_c #derivada (1-u-v+C(u,v)) em rela??o a v

#---------------------------------------------------------------------------------------------
# Likelihood function
veros=function(theta){
  
  al_t=exp(theta[1])
  ll_t=exp(theta[2])
  
  al_c=exp(theta[3])
  ll_c=exp(theta[4])
  
  beta_te=theta[5:6]
  beta_ce=theta[7:8]
  
  logl= sum(   rho*(ind * ( log(h_t(al_t,ll_t,beta_te))-H_t(al_t,ll_t,beta_te)+log(-dC_t(H_t(al_t,ll_t,beta_te),H_c(al_c,ll_c,beta_ce))+0.0000001)  )
                    +(1-ind)* ( log(h_c(al_c,ll_c,beta_ce))-H_c(al_c,ll_c,beta_ce)+log(-dC_c(H_t(al_t,ll_t,beta_te),H_c(al_c,ll_c,beta_ce))+0.0000001)))+
                 (1-rho)* log(C_a(H_t(al_t,ll_t,beta_te),H_c(al_c,ll_c,beta_ce))+0.0000001))
  
  return(-logl)
}

###---------------------------------------------------------------------------------------------###
# Fit with alpha=3 (true value)
###---------------------------------------------------------------------------------------------###

param_real=c((c(a_t,l_t,a_c,l_c)),beta_t,beta_c)   # vector with the true values of the parameters
theta <- rep(0.1,8) # initial values
alfa <- 3 # fixed alpha 
data=datasets[[1]][[2]][,1]
ind=datasets[[1]][[2]][,2]
rho=datasets[[1]][[2]][,3]
covariaveis_t=datasets[[1]][[3]]
covariaveis_c=datasets[[1]][[4]]

hessiana=matrix(NA,nrow = k,ncol = 8)
est=matrix(NA,nrow = k,ncol = 8)
loglik=matrix(NA,nrow = k,ncol = 1)

for(j in 1:k){
  
  data=datasets[[j]][[2]][,1]
  ind=datasets[[j]][[2]][,2]
  rho=datasets[[j]][[2]][,3]
  covariaveis_t=datasets[[j]][[3]]
  covariaveis_c=datasets[[j]][[4]]      
  
  tt= optim(rep(0.01, 8),veros)
  tt= optim(tt$par,veros)
  tt= optim(tt$par,veros,hessian = T)
  hessiana[j,]=sqrt(diag(solve(tt$hessian)))
  est[j,]=c(exp(tt$par[1:4]),tt$par[5:8])
  loglik[j,] = tt$value
  print(j)
}

# saving the estimated results
write.table(est,file='Est_Weibull_alpha3.txt', col.names = F, row.names = F)   #estimates saved in sequence c((c(a_t,l_t,a_c,l_c)),beta_t,beta_c)
write.table(hessiana,file='EP_Weibull_alpha3.txt', col.names = F, row.names = F)  #standard error saved in sequence c((c(a_t,l_t,a_c,l_c)),beta_t,beta_c)
write.table(loglik,file='Loglik_alpha3.txt', col.names = F, row.names = F)

###---------------------------------------------------------------------------------------------###
# Analysis of outputs, assuming alha=3
###---------------------------------------------------------------------------------------------###

est1      <- cbind(est[,1:2],est[,5:6],est[,3:4],est[,7:8])
ErroP     <- cbind(sqrt(((est[,1]^2))*hessiana[,1]^2), sqrt(((est[,2]^2))*hessiana[,2]^2), hessiana[,5:6],sqrt(((est[,3]^2))*hessiana[,3]^2),sqrt(((est[,4]^2))*hessiana[,4]^2),hessiana[,7:8]) # to alpha^T, Lambda^T, alpha^C, lambda^C delta method is done

para_real <- c(param_real[1:2],param_real[5:6],param_real[3:4],param_real[7:8])

RV <- (est1-matrix(rep(para_real,k),nrow = k,byrow =T))*100/matrix(rep(para_real,k),nrow = k,byrow =T)  
boxplot(RV)

ampitude <-  (est1+1.96*ErroP)-(est1-1.96*ErroP)    

# coverage probability
R <- nrow(est1)
lim_inf_r <- matrix(NA,nrow(est1),8)
lim_sup_r <- matrix(NA,nrow(est1),8)

for (i in 1:8){
  lim_inf_r[,i] <- est1[,i] -1.96*sd(est1[,i])
  lim_sup_r[,i] <- est1[,i] +1.96*sd(est1[,i])
}
Prob_cobertura <- NULL
for (i in 1:8){
  Prob_cobertura[i] <- (sum(lim_sup_r[,i] >= para_real[i] & lim_inf_r[,i] <= para_real[i]))/R
}

# Table
stat<-matrix(rep(NA,7*8),ncol=7,nrow=8)
stat<-as.data.frame(stat)
row.names(stat) <- c("$\\alpha^{T}$","$\\lambda^{T}$","$\\beta^{T}_{1}$","$\\beta^{T}_{2}$","$\\alpha^{C}$","$\\lambda^{C}$", "$\\beta^{C}_{1}$","$\\beta^{C}_{2}$")  #$\\alpha$: fixed
colnames(stat)<- c( "Real", "Estimate", "SE", "SD", "RV%", "Amplitude", "PC")   #SE:standard error; RV\% vicio relativo*100; Amplitude:amplitude do IC
stat[,1] <- para_real
stat[,2] <- paste(round(apply(est1[,],2,mean),3))
stat[,3] <- paste(round(apply(ErroP[,],2,mean),3)) # average of 500 standard errors
stat[,4] <- paste(round(apply(est1[,],2,sd),3))  # standard deviation
stat[,5] <- paste(round(apply(RV[,],2,mean),3)) # RV\% relative bias
stat[,6] <- paste(round(apply(ampitude[,],2,mean),3)) # Amplitude
stat[,7] <- round(Prob_cobertura,3) # coverage probability

write.table(stat,file='Tab_Weibull_alpha3.txt') 

