# Função para simular todos os métodos ao mesmo tempo
library(tidyverse)


rm(list=ls())
# Funcoes da Distribuição Maxwell
source("funcoes.R")
source("methods.R")

# SETs
set.seed(124)
nrep <- 10
mu_values <- c(0.2, 0.5, 1, 1.5)  # Valores de mu diferentes
n <- c(30, 80, 300, 600)
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

#==========/==========/==========/==========/==========/==========/==========/==========/

# Simulação com todos os metodos
simu_methods <- function(nrep, n, mu_true, methods) {
  
  mu_estim <- data.frame(matrix(ncol = length(methods), nrow = nrep))
  
    for (i in 1:nrep) {
      
      sample <- rmax(n, mu_true)
  
      for (k in 1:length(methods)) {
        meth <- methods[k]
        
        mu_estim[i, k] <- estimate_mu(sample, meth)
    }
  }
  
  colnames(mu_estim) <- methods
  
  return(mu_estim)
}


run_simulation <- function(nrep, n_values, mu_values, methods) {
  
  results <- data.frame(
    n = integer(),
    mu_true = numeric(),
    method = character(),
    mu_est = numeric()
  )
  
  # Loops para tamanhos de amostra, valores de mu e métodos
  for (n in n_values) {
    for (mu in mu_values) {
      
      mu_estimations <- simu_methods(nrep, n, mu, methods)
      
      # Armazenar os resultados
      for (meth in methods) {
        results <- rbind(results, data.frame(
          n = rep(n, nrep),
          mu_true = rep(mu, nrep),
          method = rep(meth, nrep),
          mu_est = mu_estimations[[meth]]
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

