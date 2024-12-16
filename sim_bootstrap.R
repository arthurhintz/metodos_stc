library(tidyverse)

rm(list=ls())
# Funcoes da Distribuição Maxwell
source("funcoes.R")
source("methods.R")

set.seed(124)
nrep <- 5
mu_values <- 1  # Valor de mu (pode ser expandido para mais valores)
n <- c(30, 80, 300, 600)
n_scen <- length(n)
method <- c("MM", "MLE", "AD", "MPS")  # Métodos avaliados
B <- 10  # Reamostragem de Bootstrap

# Função genérica para estimar mu usando diferentes métodos
estimate_mu <- function(x, method) {
  fn <- switch(method,
               "MM" = function(mu) mm_stat(x, mu),
               "MLE" = function(mu) loglik_max(x, mu),
               "AD" = function(mu) ad_stat(x, mu),
               "MPS" = function(mu) mps_stat(x, mu),
               stop("Método inválido"))
  
  # Estimativa do parâmetro mu usando o método especificado
  result <- optim(
    par = 1,  # Chute inicial para otimização
    fn = fn,
    method = "BFGS"
  )
  
  return(abs(result$par))  # Retorna o valor estimado de mu
}


# Função para realizar o Bootstrap e corrigir o viés de mu
bootstrap_correction <- function(amostra_original, method, B) {
  
  # Vetores para armazenar os estimadores de mu em cada iteração do Bootstrap
  estimadores_mu <- numeric(B)
  
  # Reamostragem com Bootstrap
  for (b in 1:B) {
    # Gera uma amostra bootstrap com reposição
    amostra_boot <- sample(amostra_original, size = length(amostra_original), replace = TRUE)
    
    # Estima o parâmetro mu para a amostra bootstrap
    estimadores_mu[b] <- estimate_mu(amostra_boot, method)
  }
  
  # Calcula a média dos estimadores bootstrap
  media_mu <- mean(estimadores_mu)
  
  # Estima o viés (diferença entre a média bootstrap e a estimativa original)
  viés_mu <- media_mu - mean(amostra_original)
  
  return(list(
    media_mu = media_mu,
    viés_mu = viés_mu
  ))
}

# Função principal para rodar as simulações Monte Carlo
run_simulations <- function(nrep, n_values, mu_true, methods, B) {
  results <- data.frame()  # Inicializa um data frame para armazenar os resultados
  
  for (n in n_values) {
    for (method in methods) {
      
      # Vetores para armazenar as médias e viéses de cada repetição
      media_mc_mu <- numeric(nrep)
      viés_mc_mu <- numeric(nrep)
      
      for (mc in 1:nrep) {
        # Gerar a amostra original usando a função rmax
        amostra_original <- rmax(n, mu_true) 
        
        cat('n = ', n, 'metodo = ', method, ' mc =', mc)
        
        # Realizar o Bootstrap para corrigir o viés
        resultados_boot <- bootstrap_correction(amostra_original, method, B)
        
        # Armazenar os resultados de cada Monte Carlo
        media_mc_mu[mc] <- resultados_boot$media_mu
        viés_mc_mu[mc] <- resultados_boot$viés_mu
      }
      
      # Adiciona os resultados ao data frame
      results <- rbind(results, data.frame(
        n = rep(n, nrep),
        mu_true = rep(mu_true, nrep),
        method = rep(method, nrep),
        media_mu = media_mc_mu,
        viés_mu = viés_mc_mu
      ))
    }
  }
  
  return(results)
}

set.seed(124)
# Rodar as simulações

nrep <- 500
mu_values <- 1  # Valor de mu (pode ser expandido para mais valores)
n <- c(30, 80, 250, 500)
n_scen <- length(n)
method <- c("MM", "MLE", "AD", "MPS")  # Métodos avaliados
B <- 5000 # Reamostragem de Bootstrap

resultados <- run_simulations(nrep, n, mu_values, method, B)
resultados

head(resultados[resultados$method == "MLE" & resultados$n == 30, ])  

write.csv(resultados, "resultados_simulacao2.csv", row.names = FALSE)


