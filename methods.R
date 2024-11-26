source("funcoes.R")

# Métodos dos Momentos
mm_stat <- function(x, mu) {
  media = mean(x)

  return((media - mu)^2)
}

# MV
loglik_max <- function(x, mu) {
  
  mu <- mu[1]  
  
  log_lik <- sum(log(dmax(x, mu)))
  
  return(-log_lik)
}


# Anderson-Darling

ad_stat <- function(x, mu) {
  n <- length(x)
  sorted_data <- sort(x)  # Ordenar os dados
  cdf_values <- pmax(sorted_data, mu) 
  
  # Fórmula da estatística Anderson-Darling
  A2 <- -n - (1 / n) * sum(
    (2 * (1:n) - 1) * (log(cdf_values) + log(1 - rev(cdf_values)))
  )
  return(A2)  # Retorna o valor da estatística A^2
}


# Log do produto dos espaçamentos (MPS)

mps_stat <- function(x, mu) {
  n <- length(x)
  sorted_data <- sort(x)  # Ordenar os dados
  cdf_values <- pmax(sorted_data, mu)  
  
  # Espaçamentos (incluir limites 0 e 1)
  spacings <- c(cdf_values, 1) - c(0, cdf_values)
  spacings[spacings <= 0] <- 1e-10  # Evitar log(0)
  
  # Retorna o log negativo do produto dos espaçamentos
  return(-sum(log(spacings)))
}


