# Função para simular todos os métodos ao mesmo tempo
library(tidyverse)


rm(list=ls())
# Funcoes da Distribuição Maxwell
source("funcoes.R")
source("methods.R")

# SETs
set.seed(1234)
nrep <- 5
mu_values <- c(0.2, 0.5, 1, 1.5)  # Valores de mu diferentes
n <- c(30, 100, 300, 600)
n_scen <- length(n)
method <- c("MM","MLE", "AD", "MPS")  # Métodos avaliados


# Função genérica para estimar mu usando diferentes métodos
estimate_mu <- function(x, method) {
  fn <- switch(method,
               "MM" = function(mu) mm_stat(x, mu),
               "MLE" = function(mu) loglik_max(x, mu),
               "AD" = function(mu) ad_stat(x, mu),
               "MPS" = function(mu) mps_stat(x, mu),
               stop("Método inválido"))
  result <- optim(
    par = 1,  # Chute inicial
    fn = fn,
    method = "BFGS"
  )
  return(result$par)  # Retorna o valor estimado de lambda
}


# Função para simular estimativas x
simu_methods <- function(nrep, n, mu_true, method) {
  
  mu_estim <- numeric(nrep)
  
  for (i in 1:nrep) {
    sample <- rmax(n, mu_true)
    mu_estim[i] <- estimate_mu(sample, method)
  }
  return(mu_estim)
}


#==========/==========/==========/==========/==========/==========/==========/==========/

# Função principal de simulação
run_simulation <- function(nrep, n_values, mu_values, methods) {

    # DataFrame para armazenar resultados
  results <- data.frame(
    n = integer(),
    mu_true = numeric(),
    method = character(),
    mu_est = numeric()
  )
  
  # Loops para tamanhos de amostra, valores de mu e métodos
  for (n in n_values) {
    for (mu in mu_values) {
      for (meth in methods) {

        mu_estim <- simu_methods(nrep, n, mu, meth)
        
        # Armazenar os resultados 
        results <- rbind(results, data.frame(
          n = rep(n, nrep),
          mu_true = rep(mu, nrep),
          method = rep(meth, nrep),
          mu_est = mu_estim
        ))
      }
    }
  }
  
  return(results)
}


results <- run_simulation(nrep, n, mu_values, method)

write.table(
  results,
  file = "simulation_results.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# Análise gráfica
ggplot(results, aes(x = factor(n), y = mu_est, fill = method)) +
  geom_boxplot() +
  facet_wrap(~ mu_true, scales = "free_y") +
  labs(
    title = "Estimativa de mu para diferentes tamanhos de amostra e métodos",
    x = "Tamanho da amostra",
    y = "Estimativa de mu"
  ) +
  theme_minimal()

# Análise pontual
 
resu <- results |> 
  group_by(n, mu_true, method) |> 
  summarise(media = mean(mu_est))

