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
  ext <- uniroot(function(x) MxARMA::pmax(x, mu) - u, interval = c(0, 1), tol = 1e-300, extendInt = "yes")$root
  return(abs(ext))
}
#==========/==========/==========/==========/==========/==========/==========/==========/
# Função Geradora

rmax <- function(n, mu) {
  u <- runif(n)
  r <- rep(NA, n)
  for (i in 1:n) {
    r[i] <- MxARMA::qmax(u[i], mu)
  }
  return(r)
}
#==========/==========/==========/==========/==========/==========/==========/==========/
# Log Verossimilhança

# estão diferentes as funções conferir isso

log_vero_max <- function(y, mu){

  ll <- log(32) - 2*log(pi) - 3*log(mu) + log(y^2) 
        + ((-4*y^2)/(pi*mu^2))
  
  ll <- -sum(ll)
  
  return(ll)
}
#ou

loglik_max <- function(y, mu) {
  
  mu <- mu[1]  
  
  log_lik <- sum(log(dmax(y, mu)))
  
  return(-log_lik)
}

#==========/==========/==========/==========/==========/==========/==========/==========/
set.seed(1248)
mu_true <- 2    
n <- 100      # Tamanho da amostra

x_simulado <- rmax(n = n, mu = mu_true)

# 5. Estimar os Parâmetros via Máxima Verossimilhança

valores_iniciais <- 1

# Usamos a função 'optim' para encontrar os valores que maximizam a log-verossimilhança
resultado <- optim(par = valores_iniciais, 
                   fn = loglik_max, y = x_simulado,
                   method = "SANN")

resultado$par  

loglik_max(y = 2, mu = 2)
log_vero_max(y = 2, mu = 2)


