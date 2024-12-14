source("funcoes.R")


# Tamanho da imagem
{width <- 3
  height <- 3
  mar_b<-2.5
  mar_e<-2.5
  mar_c<-0.5
  mar_d<-0.5
  dist_text<-1.5
  dist_tick<-0.5
}


library(ggplot2)


fromx <- 0
tox <- 2.5
yliminf <- 0
ylimsup <- 5

# Valores de mu
mu_values <- c(0.2, 0.5, 0.8, 1.2)
comp <- length(mu_values)


# Cores e tipos de linha
colors <- 1:comp
line_types <- 1:comp


# Create a data frame with x values and corresponding densities for each mu
x_values <- seq(fromx, tox, length.out = 1000) # Generate x values
density_data <- data.frame(
  x = rep(x_values, times = length(mu_values)),
  y = unlist(lapply(mu_values, function(mu) dmax(x_values, mu))),
  mu = rep(mu_values, each = length(x_values))
)

# Convert mu to a factor for color and linetype mapping
density_data$mu <- factor(density_data$mu, levels = mu_values)

# Create the plot
pdf(file = "fdpMW.pdf",width = width, height = height,family = "Times")
par(mar=c(mar_b, mar_e, mar_c, mar_d)) 
par(mgp=c(dist_text, dist_tick, 0))
ggplot(density_data, aes(x = x, y = y, color = mu, linetype = mu)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = line_types) +
  labs(y = expression("f(y)"), x = "x") +
  theme_minimal(base_size = 15) +
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2)),
         linetype = guide_legend(override.aes = list(size = 2)))
dev.off()



# Define the dmax function (ensure this is defined in your environment)
dmax <- function(x, mu) {
  # Example implementation, replace with your actual function
  exp(-x / mu) * (x >= 0) / mu
}

# Generate x-values for the plot
x_values <- seq(fromx, tox, length.out = 1000)

# Create a data frame for the plot
density_data <- data.frame(
  x = rep(x_values, times = length(mu_values)),
  y = unlist(lapply(mu_values, function(mu) dmax(x_values, mu))),
  mu = rep(mu_values, each = length(x_values))
)

# Convert mu to a factor for use in ggplot aesthetics
density_data$mu <- factor(density_data$mu, levels = mu_values)

pdf(file = "fdaMW.pdf",width = width, height = height,family = "Times")
par(mar=c(mar_b, mar_e, mar_c, mar_d)) 
par(mgp=c(dist_text, dist_tick, 0))
ggplot(density_data, aes(x = x, y = y, color = mu, linetype = mu)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = line_types) +
  labs(y = expression("f(y)"), x = "x") +
  theme_minimal(base_size = 15) +
  coord_cartesian(ylim = c(yliminf, ylimsup)) + # Set y-axis limits
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2)),
         linetype = guide_legend(override.aes = list(size = 2)))
dev.off()



fromx <- 0
tox <- 2.5
yliminf <- 0
ylimsup <- 5

# Valores de mu
mu_values <- c(0.2, 0.5, 0.8, 1.2)
comp <- length(mu_values)


# Cores e tipos de linha
colors <- 1:comp
line_types <- 1:comp

curve(dmax(x, mu = mu_values[1]), from = fromx, to = tox, 
        lty = line_types[1], type = "l", 
        ylab = expression("f(y)"), ylim = c(yliminf, ylimsup), 
        col = colors[1], lwd = 2.0, cex.lab = 1.5)
  
  for (i in 2:length(mu_values)) {
    mu <- mu_values[i]
    curve(dmax(x, mu), from = fromx, to = tox, add = TRUE, lty = line_types[i], 
          col = colors[i], lwd = 2.0)
  }
  
  legend(x = (tox - 0.5), y = ylimsup, legend = sapply(mu_values, function(x) bquote(mu == .(x))),
         col = colors, lty = line_types, lwd = rep(4, comp), bty = "n", cex = 1.2)

  
  # Define the function to plot each curve
  plot_curve <- function(mu, lty, col) {
    curve(dmax(x, mu), from = fromx, to = tox, add = TRUE, lty = lty, col = col, lwd = 2.0)
  }


fromx <- 0
tox <- 2.5
yliminf <- 0
ylimsup <- 1

# Valores de mu
mu_values <- c(0.2, 0.5, 0.8, 1.2)
comp <- length(mu_values)


# Cores e tipos de linha
colors <- 1:comp
line_types <- 1:comp

{
  
  curve(pmax(x, mu = mu_values[1]), from = fromx, to = tox, 
        lty = line_types[1], type = "l", 
        ylab = expression("f(y)"), ylim = c(yliminf, ylimsup), 
        col = colors[1], lwd = 2.0, cex.lab = 1.5)
  
  for (i in 2:length(mu_values)) {
    mu <- mu_values[i]
    curve(pmax(x, mu), from = fromx, to = tox, add = TRUE, lty = line_types[i], 
          col = colors[i], lwd = 2.0)
  }
  
  legend(x = (tox - 0.5), y = ylimsup, legend = sapply(mu_values, function(x) bquote(mu == .(x))),
         col = colors, lty = line_types, lwd = rep(4, comp), bty = "n", cex = 1.2)
}


fromx <- 0
tox <- 10
yliminf <- -100
ylimsup <- 1

# Valores de mu
x_values <- c(1,4,8,12)
comp <- length(x_values)


# Cores e tipos de linha
colors <- 1:comp
line_types <- 1:comp

mu_values <- seq(fromx, tox, length.out = 100)

y_values <- sapply(x_values, function(x) sapply(mu_values, loglik_max, x = x))

# Plotagem
plot(mu_values, y_values[, 1], type = "l", lty = line_types[1], col = colors[1], 
     lwd = 2.0, xlab = expression(mu), ylab = expression("l(mu;y)"), 
     ylim = c(yliminf, ylimsup), cex.lab = 1.2)

# Adicionar as outras curvas
for (i in 2:comp) {
  lines(mu_values, y_values[, i], lty = line_types[i], col = colors[i], lwd = 2.0)
}

# Adicionar legenda
legend(x = (tox - 2.5), y = (ylimsup - 10), legend = sapply(x_values, function(x) bquote(x == .(x))),
       col = colors, lty = line_types, lwd = rep(4, comp), bty = "n", cex = 0.8)