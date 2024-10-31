# Simulação de Monte Carlo

rm(list=ls(all=TRUE))
# Funcoes
source("funcoes.R")
# SET
set.seed(1234)
nrep <- 5000
mu_values <- c(0.2, 0.5, 1, 1.5, 2)  # Valores de mu diferentes
chute <- 1
n <- seq(50, 1000, 50)
n_scen <- length(n)
alpha <- 0.05
medidas <- c("mean", "RB(%)", "MSE", "CR")

resultados_mu <- list()

for (mu in mu_values) {
  
  LI <- qmax(alpha, mu)
  LS <- qmax((1 - alpha), mu)
  
  resu <- matrix(NA, ncol = n_scen, nrow = nrep)
  colnames(resu) <- c(as.character(n))
  
  # SIMU
  for (k in seq_along(n)){
    
    for (i in 1:nrep){
      
      x_sim <- rmax(n = k, mu = mu)
      
      resu[i, k] <- optim(par = chute, 
                          fn = loglik_max, x = x_sim,
                          method = "SANN")$par
    }
  }
  
  par_est <- apply(resu, 2, mean)
  vviesm <- (par_est - mu)
  viest <- (par_est - mu) / mu * 100
  vvarm <- apply(resu, 2, var)
  veqmm <- vviesm^2 + vvarm
  txcobertura <- colMeans(resu >= LI & resu <= LS)
  
  resu_final <- cbind(round(par_est, digits = 3),
                      round(viest, digits = 3), 
                      round(veqmm, digits = 3),
                      round(txcobertura, digits = 3))
  
  colnames(resu_final) <- medidas
  
  resultados_mu[[paste0("mu_", mu)]] <- resu_final
}


for (mu_name in names(resultados_mu)) {
  write.table(resultados_mu[[mu_name]], 
              paste0(mu_name, "_resultado.txt"))
}

