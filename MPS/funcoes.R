# Maxwell Reparametrizada
# Função Desnsidade de probabilidade

dmax <- function(y, mu) {
  ((32 * y^2) / (pi^2 * mu^3)) * exp(-4 * y^2 / (pi * mu^2))
}
#==========/==========/==========/==========/==========/==========/==========/==========/

# Funcão Acumulada
pmax <- function(y, mu) {
  (2/sqrt(pi)) * gamma(3/2) * pgamma((4 * y^2)/(pi * mu^2), 3/2, lower.tail = TRUE)
}
#==========/==========/==========/==========/==========/==========/==========/==========/
# Função quantílica

qmax <- function(u, mu) {
  ext <- uniroot(function(x) pmax(x, mu) - u, interval = c(0, 1), tol = 1e-300, extendInt = "yes")$root
  return(abs(ext))
}
#==========/==========/==========/==========/==========/==========/==========/==========/
# Função Geradora de números aleatórios

rmax <- function(n, mu) {
  u <- runif(n)
  r <- rep(NA, n)
  for (i in 1:n) {
    r[i] <- qmax(u[i], mu)
  }
  return(r)
}
#==========/==========/==========/==========/==========/==========/==========/==========/
# Log Verossimilhança MLE

# estão diferentes as funções conferir isso

log_vero_max <- function(x, mu){
  mu <- mu[1]  
  
  ll <- log(32) + 2*log(x) - (2*log(pi) + 3*log(mu)) +  
        + ((-4*x^2)/(pi*mu^2))
  
  ll <- sum(ll)
  
  return(-ll)
}
#ou

loglik_max <- function(x, mu) {
  
  mu <- mu[1]  
  
  log_lik <- sum(log(dmax(x, mu)))
  
  return(-log_lik)
}

#==========/==========/==========/==========/==========/==========/==========/==========/
# Vetor Escore MLE

Score_vector <- function(x){
  
  mu_hat <- sqrt(8*sum(x**2)/(3*pi*length(x)))
  
  return(mu_hat)
}

## MPS estimator functions -----------------------------------------------------

MPS_Mw <- function(teta,y,m.optim = 1,Geo_mean = T)
{
  mu <- teta[1]

  y <- sort(y)
  tmp <-pmax(y = y, mu = mu)
  Fi <- c(0,tmp)
  Ff <- c(tmp,1)
  D <- Ff-Fi
  for(j in 1:length(D)){
    if(is.null(D[j])){
      D[j]<-dmax(y = y, mu = mu)
    }
  }
  if(Geo_mean == T){
    H<-suppressWarnings(mean(log(D)))
    if(m.optim==-1){return(-H)} 
    if(m.optim==1){return(H)}
  }else{
    if(m.optim==-1){return(-D)} 
    if(m.optim==1){return(D)}}
}

## Log-likelihood type MPS function --------------------------------------------

# llike_MPS_Mw <- function(x,y)
# {
#   log_like <- sum(
#     log(x[2]*x[1]*x[3])-2*log(exp((1-exp(-log(y))^x[2])*x[1])*(x[3]-1)+1)+
#       (1-exp(-log(y))^x[2])*x[1]+log(-log(y))*(x[2]-1)+(-log(y))^x[2]-log(y)
#   )
#   return(log_like)
# }

# Score vector type MPS function -----------------------------------------------

vscore_MPS_Mw <- function(teta,y,m.optim=1,Geo_mean=F)
{
  D <- MPS_Mw(teta=teta,y=y,m.optim = m.optim,Geo_mean=F)
  mu <- teta[1]
  n <- length(y)
  y <- sort(y)
  
  # derivada em relação a mu
  wt <- - (32*(y**3)*exp(4*(y**2)/(pi*(mu**2)))) / ((pi**2)*(mu**2))
  
  TFi <- matrix(c(c(0,wt),c(0,rt),c(0,st)),nrow = 3, ncol = n+1,byrow = T)
  TFf <- matrix(c(c(wt,0),c(rt,0),c(st,0)),nrow = 3, ncol = n+1,byrow = T)
  TD <- TFf-TFi
  for(j in 1:(length(TD)/3)){
    if(is.null(TD[1,j])){
      TD[1,j] <- sort(dmax(y = y, mu = mu))[j]
    }
  }
  V <- suppressWarnings(colMeans(t(TD)/D))
  
  if(m.optim==-1){return(-V)} 
  if(m.optim==1){return(V)} 
} 

#### ---------------------------------------------------------------------- ####