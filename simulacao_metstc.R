# Simulação de Monte Carlo

rm(list=ls(all=TRUE))
# Funcoes
source("funcoes.R")

# SET
set.seed(1234)
nrep <- 5
mu <- 2
chute <- 1
n <- seq(50, 1000, 50)
n_scen <- length(n)
alpha <- 0.05

LI <- qmax(alpha, mu)
LS <- qmax((1 - alpha), mu)


medidas <- c("mean", "RB(%)", "MSE", "CR")

resu <- matrix(NA, ncol = n_scen, nrow = nrep)
colnames(resu) <- c(as.character(n))

# SIMU
for (k in seq_along(n)){
  
  for (i in 1:nrep){
    
    x_sim <- rmax(n = k, mu = mu)
    
    resu[i,k] <- optim(par = chute, 
                       fn = loglik_max, x = x_sim,
                       method = "SANN")$par
  }
}
resu

par_est <- apply(resu, 2, mean)
vviesm <- (par_est - mu)
viest <- (par_est - mu) / mu * 100
vvarm <- apply(resu, 2, var)
veqmm <- vviesm^2 + vvarm
txcobertura <- colMeans(resu >= LI & resu <= LS)



resu <- cbind(round(par_est, digits = 3),
              round(viest, digits = 3), 
              round(veqmm, digits = 3),
              round(txcobertura, digits = 3))


colnames(resu) <- medidas

write.table(resu, "resultado.txt")





