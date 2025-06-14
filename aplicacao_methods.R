source("methods.R")
source("funcoes.R")

#url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv"
#wine <- read.csv(url, sep = ";")



dados <- readxl::read_xlsx("base - Copia.xlsx")


k <- ncol(dados)

summary(dados)

dado <- dados$F

dado <- na.omit(dado)

hist(dado)

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

samples <- c(length(dado))

results <- array(NA, dim = c(length(methods), 4, length(samples)),
                 dimnames = list(methods, c("Estimativa", "Media_Original", "Vies", "EQM"), as.character(samples)))

# Loop corrigido
for (i in seq_along(samples)) {
  n <- samples[i]
  x <- sample(dado, n)
  mean_x <- mean(x)
  
  for (j in seq_along(methods)) {
    m <- methods[j]
    resu <- estimate_mu(x, method = m)
    mu_hat <- abs(resu$par)
    
    if (is.matrix(resu$hessian) && det(resu$hessian) != 0) {
      vcov <- solve(resu$hessian)
      ep <- sqrt(diag(vcov))
    } else {
      ep <- NA
    }
    
    bias <- mu_hat - mean_x
    eqm <- ep^2 + bias^2
    
    results[j, , i] <- c(mu_hat, mean_x, bias, eqm)
  }
}


print(results)

#==========/==========

# Cores nomeadas
cores <- c("blue", "red", "orange", "purple")
names(cores) <- methods

xmin <- min(dado)
xmax <- max(dado)
xlim <- seq(xmin, xmax, 0.01)

# Histograma com densidade empírica
hist(x, breaks = 5, probability = T, col = "lightgray",
     main = "Comparação das curvas estimadas", xlab = "x")

# Adiciona curvas teóricas dmax para cada estimativa de mu
for (m in methods) {
  mu_est <- results[m, "Estimativa", 1]
  
  curve(dmax(x, mu = mu_est), from = min(x), to = max(x),
        col = cores[m], lwd = 2, add = TRUE)
}

# Adiciona densidade empírica
lines(density(x), col = "darkgreen", lwd = 2, lty = 2)

# Legenda
legend("topright",
       legend = c(paste("dmax -", methods), "Densidade"),
       col = c(cores, "darkgreen"),
       lty = c(rep(1, length(methods)), 2),
       lwd = 2,
       cex = 0.5)  


# Usar artigo para me basear
#usar bootstrap
# Usar dois tipo de amostragem:
# ranked set sampling (RSS) and simple random sampling (SRS) designs. 
# Fazer o grafico das probabilidades
# faz sentido incluir um método bayesiano ?





