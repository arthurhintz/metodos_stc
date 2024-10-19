# Maxwell Reparametrizada
# Função Desnsidade de probabilidade

dmax <- function(y, mu) {
  ((32 * y^2) / (pi^2 * mu^3)) * exp(-4 * y^2 / (pi * mu^2))
}

# Funcão Acumulada
pmax <- function(y, mu) {
  (2/sqrt(pi)) * gamma(3/2) * pgamma((4 * y^2)/(pi * mu^2), 3/2, lower.tail = TRUE)
}

# Função quantílica

qmax <- function(u, mu) {
  ext <- uniroot(function(x) MxARMA::pmax(x, mu) - u, interval = c(0, 1), tol = 1e-300, extendInt = "yes")$root
  return(abs(ext))
}

# Função Geradora
rmax <- function(n, mu) {
  u <- runif(n)
  r <- rep(NA, n)
  for (i in 1:n) {
    r[i] <- MxARMA::qmax(u[i], mu)
  }
  return(r)
}

# Log Verossimilhança

log_vero_max <- function(y, mu){

  ll <- log(32) - 2*log(pi) - 3*log(mu) + log(y^2) 
        + ((-4*y^2)/(pi*mu^2))
  
  sum(ll)
}


