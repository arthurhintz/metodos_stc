source("methods.R")
source("funcoes.R")

url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv"
wine <- read.csv(url, sep = ";")
head(wine)

x <- wine$pH

x <- sample(x, 70)

hist(x)

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
    method = "SANN",
    hessian = T
  )
  return(result) 
}


methods <- c("MM", "MLE", "AD", "MPS")

results <- list()

n <- length(x)

mean_x <- mean(x)

for (m in methods) {
  resu <- estimate_mu(x, method = m)
  
  mu_hat <- resu$par
  
  vcov <- solve(resu$hessian)
  
  ep <- sqrt(diag(vcov))
  
  
  bias <- mu_hat - mean_x
  eqm <- ep^2 + bias
  
  results[[m]] <- c(
    Estimativa = mu_hat,
    Média_Original = mean_x,
    Viés = bias,
    EQM = eqm
  )
}

result_df <- do.call(rbind, results)
result_df <- as.data.frame(result_df)

result_df[] <- lapply(result_df, function(col) as.numeric(as.character(col)))

result_df <- round(result_df, 4)
print(result_df)

#==========/==========

# Cores nomeadas
cores <- c("blue", "red", "orange", "purple")
names(cores) <- methods

# Histograma com densidade empírica
hist(x, breaks = 10, probability = TRUE, col = "lightgray",
     main = "Comparação das curvas estimadas", xlab = "pH")

# Adiciona curvas teóricas dmax para cada estimativa de mu
for (m in methods) {
  mu_est <- results[[m]][["Estimativa"]]
  
  curve(dmax(x, mu = mu_est), from = min(x), to = max(x),
        col = cores[m], lwd = 2, add = TRUE)
}

# Adiciona densidade empírica
lines(density(x), col = "darkgreen", lwd = 2, lty = 2)

# Legenda
legend("topright",
       legend = c(paste("dmax -", methods), "Densidade empírica"),
       col = c(cores, "darkgreen"),
       lty = c(rep(1, length(methods)), 2),
       lwd = 2)

# Usar artigo para me basear
#usar bootstrap
# Usar dois tipo de amostragem:
# ranked set sampling (RSS) and simple random sampling (SRS) designs. 
# Fazer o grafico das probabilidades
# faz sentido incluir um método bayesiano ?





