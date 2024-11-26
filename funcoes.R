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
# set.seed(1248)
# mu_true <- 2    
# n <- 1000      # Tamanho da amostra
# 
# x_simulado <- rmax(n = n, mu = mu_true)
# 
# # 5. Estimar os Parâmetros via Máxima Verossimilhança
# 
# valores_iniciais <- 5
# 
# # Usamos a função 'optim' para encontrar os valores que maximizam a log-verossimilhança
# resultado <- optim(par = valores_iniciais, 
#                    fn = loglik_max, x = x_simulado,
#                    method = "SANN")
# 
# #2
# resultado <- optim(par = valores_iniciais, 
#                    fn = log_vero_max, x = x_simulado,
#                    method = "SANN")
# 
# 
#==========/==========/==========/==========/==========/==========/==========/==========/
# Informação de Fisher

fisher <- function(n, mu){
  
  a <- n * 6 / mu^2 
  
  return(a)
  
}

