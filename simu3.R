# simu 3


rm(list=ls(all=TRUE))
# Funcoes
source("funcoes.R")
# SET
set.seed(1234)
nrep <- 10
mu_values <- c(0.2, 0.5, 1, 1.5)  # Valores de mu diferentes
chute <- 1
n <- c(50, 150, 300, 600)
n_scen <- length(n)
z <- qnorm(0.975)
medidas <- c("mean", "RB(%)", "MSE", "CR")

resultados_mu <- list()

for (mu in mu_values) {
  
  resu <- matrix(NA, ncol = n_scen, nrow = nrep)
  colnames(resu) <- c(as.character(n))
  
  ic <- matrix(0, nrep, n_scen)
  # SIMU
  for (k in 1:n_scen){
    
    p <- n[k]
    
    for (i in 1:nrep){
      
      x <- rmax(n = p, mu = mu)
      
      muhat <- sqrt((8*sum(x^2)/(p*3*pi)))
      
      resu[i,k] <- muhat
      
      ifsher <- p * 6 / muhat^2
      ep <- sqrt(1/ifsher)
      
      LI <- muhat - z * ep
      LS <- muhat + z * ep
      
      ic[i, k] <- as.numeric((mu > LI) & (mu < LS))
      
    }
  }
  
  par_est <- apply(resu, 2, mean)
  vviesm <- (par_est - mu)
  viest <- (par_est - mu) / mu * 100
  vvarm <- apply(resu, 2, var)
  veqmm <- vviesm^2 + vvarm
  txcobertura <- apply(ic, 2, mean)
  
  resu_final <- cbind(round(par_est, digits = 3),
                      round(viest, digits = 3), 
                      round(veqmm, digits = 3),
                      round(txcobertura, digits = 3))
  
  colnames(resu_final) <- medidas
  
  resultados_mu[[paste0("mu_", mu)]] <- resu_final
}


for (mu_name in names(resultados_mu)) {
  write.table(resultados_mu[[mu_name]], 
              paste0(mu_name, "_resultado2.txt"))
}

