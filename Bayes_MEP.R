
###---------------------------------------------------------------------------------------------###
# Functions for the Clayton copula model with piecewise exponential marginal distribution,        #              # 
# with Bayesian inference                                                                         #
###---------------------------------------------------------------------------------------------###

rm(list=ls(all=TRUE)) 

wd <- getwd()
require(R2jags)
require(rjags)

###---------------------------------------------------------------------------------------------###
# Piecewise exponential functions
###---------------------------------------------------------------------------------------------###

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

### --------------------------------------------------------
hmean <- function(x)
{
  return(1/mean(1/x))
}

### --------------------------------------------------------
CPO <- function(loglik, n)
{  
  lik <- apply(loglik, 2, exp)
  CPO <- apply(lik, 2, hmean)
  LPML <- sum(log(CPO))
  aLPML <- mean(log(CPO))
  return(c(LPML, aLPML))  
}

### --------------------------------------------------------
DIC <- function(loglik, n)
{
  D <- apply( -2*loglik, 1, sum);
  pD <- 0.5*var(D)
  DIC <- mean(D) + pD;
  return( matrix(c(DIC, pD), ncol=2) )  
}

### --------------------------------------------------------
WAIC <- function(loglik, n)
{
  lpd <- sum( log( apply(exp(loglik), 2, mean) ) )
  pD <- sum( apply(loglik,  2, var) )
  WAIC <- lpd - pD
  return( matrix(c(WAIC, pD),ncol=2) )  
}

### --------------------------------------------------------
moda.dens <- function(dn, plotit=TRUE){
  ## dn: um objeto do uso da fun??o density()
  ini <- dn$x[which.max(dn$y)]
  fx <- approxfun(dn$x, dn$y)
  op <- optim(c(ini), fx, method="BFGS", control=list(fnscale=-1))
  if(plotit==TRUE){
    plot(dn)
    abline(v=op$par, col=2)
  }
  return(moda=op$par)
}

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
# Adjustment
###---------------------------------------------------------------------------------------------###

# chain values
n.chain <- 1
burnin  <- 1000 
lag <- 30
Npost <-  1000 
Nsim <- burnin + (Npost)*lag
mcmc <- list(n.chain=n.chain, burnin=burnin, lag=lag, Nsim=Nsim)

# setting hyperparameters: 
a.lambdaT = 0.01 
b.lambdaT = 0.01
a.lambdaC = 0.01
b.lambdaC = 0.01
#tau.alpha = 10^(-3)
b.alpha = 20

p = 2
q = 2
mu.t <- c(0,0)
mu.c <- c(0,0)
tau.beta <- 10^(-3)
Tau.t <- matrix(0,2,2)
Tau.c <- matrix(0,2,2)
Tau.t[1,1] <- tau.beta
Tau.t[2,2] <- tau.beta
Tau.c[1,1] <- tau.beta
Tau.c[2,2] <- tau.beta

const = 10000

bmax = 7


for(j in 334:k){
  
setwd(wd)
# ler os dados 
time=datasets[[j]][[2]][,1] #tempo
ind=datasets[[j]][[2]][,2]
rho=datasets[[j]][[2]][,3]
x_t=datasets[[j]][[3]]
x_c=datasets[[j]][[4]]
p = ncol(x_t)
q = ncol(x_c)
n = length(time)

delta.t <- ind
delta.c <- (1-ind)*rho

a <- time.grid(time,delta.t,bmax)
a <- a[[1]]
s <- time.grid(time,delta.c,bmax) 
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


# fit with JAGS
data <- list(rho=rho, ind=ind, n=n, dummy=rep(0,n), X_T=x_t ,X_C=x_c, p=p, q=q,b=b, c=c,
             id.falha=id.falha, id.cens=id.cens, xi.falha=xi.falha, xi.cens=xi.cens,
             a.lambdaT=a.lambdaT, b.lambdaT=b.lambdaT, a.lambdaC=a.lambdaC, b.lambdaC=b.lambdaC,
             mu.t=mu.t, mu.c=mu.c, Tau.t=Tau.t, Tau.c=Tau.c,
             b.alpha=b.alpha, const=const)
inits <- function(){list(beta.t= rep(0.1,p), beta.c=rep(0.1,q), lambda.t= rep(1,b), lambda.c=rep(1,c), alfa=0.5)}
parameters <- c("beta.t", "beta.c", "alfa", "lambda.t", "lambda.c", "loglik.total")
mod <- file.path(wd, "Clayton_MEP_cov.txt")

saida <- jags(data=data, inits=inits, parameters=parameters, n.iter=mcmc$Nsim, n.thin=mcmc$lag, n.burnin=burnin, n.chains=mcmc$n.chain, model.file=mod, working.directory=wd)    
amostra <- as.mcmc(saida)
par <- amostra[, -c( grep("loglik.total", varnames(amostra)), grep("deviance", varnames(amostra)))]   
par.hat <- matrix(summary(par)$statistics[,"Mean"],nrow=1)  
mediana <- matrix(summary(par)$quantiles[,"50%"],nrow=1)  
sd.par <- matrix(summary(par)$statistics[,"SD"],nrow=1)
hpd <- matrix(unlist(HPDinterval(par)), nrow=1)  

parametros <- matrix(0,ncol= length(par.hat) , nrow=length(unlist(par[,1])))
moda_para <- rep(0,length(par.hat))
for (i in 1:length(par.hat)){
  parametros[,i] <- (unlist(par[,i]))  
  dn <- density((parametros[,i]), kernel="triangular", bw=0.01)
  moda <- moda.dens(dn, plotit=FALSE)
  moda_para[i] <- moda 
}  
moda_param_post <- matrix(moda_para, nrow=1) 

par_vicio <- amostra[, -c( grep("loglik.total", varnames(amostra)), grep("deviance", varnames(amostra)),grep("lambda.c", varnames(amostra)),grep("lambda.t", varnames(amostra)))]   
value_parameters <- c( 3, beta_c, beta_t) #alfa  beta.c[1] beta.c[2]  beta.t[1] beta.t[2]  scale.c   scale.t  shape.c  shape.t
vicio.relat <- matrix(((summary(par_vicio)$statistics[,"Mean"]- value_parameters)/abs(value_parameters))*100,nrow=1)  #vicio relativo*100
nef <- matrix(effectiveSize(par), nrow=1)
loglik <- amostra[, paste("loglik.total[",1:n,"]", sep="")]
loglik <- - as.matrix(loglik, ncol=n, byrow=TRUE) + const
CPO.out <- matrix(CPO(loglik, n), nrow=1) #(LPML, aLPML)
DIC.out <- DIC(loglik, n)
WAIC.out <- WAIC(loglik, n)

wd.results_cov <- paste(wd, "/results_7int", sep="")
setwd(wd.results_cov)
write.table(par.hat, file ="media.txt", row.names=FALSE, col.names=FALSE, append=TRUE);
write.table(mediana, file ="mediana.txt", row.names=FALSE, col.names=FALSE, append=TRUE);
write.table(sd.par, file ="sd.txt", row.names=FALSE, col.names=FALSE, append=TRUE);
write.table(hpd, file ="hpd.txt", row.names=FALSE, col.names=FALSE, append=TRUE);
write.table(moda_param_post, file ="moda.txt", row.names=FALSE, col.names=FALSE, append=TRUE);
write.table(vicio.relat, file ="bias.txt", row.names=FALSE, col.names=FALSE, append=TRUE);

write.table(nef, file ="nef.txt", row.names=FALSE, col.names=FALSE, append=TRUE);
write.table(CPO.out, file ="CPO.txt", row.names=FALSE, col.names=FALSE, append=TRUE);
write.table(DIC.out, file ="DIC.txt", row.names=FALSE, col.names=FALSE, append=TRUE);
write.table(WAIC.out, file ="WAIC.txt", row.names=FALSE, col.names=FALSE, append=TRUE);


}








