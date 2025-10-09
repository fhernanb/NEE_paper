library(gamlss)
library(RelDists)

# The parameters ----------------------------------------------------------
# Combination of section 4.2 Simulation Study
true_mu    <- 2.5
true_sigma <- 25

# Useful functions to the simulation study --------------------------------

# Funcion para obtener mu_hat y sigma_hat para un valor fijo de n
simul_one <- function(size) {
  y <- rNEE(n=size, mu=true_mu, sigma=true_sigma)
  mod <- NULL
  mod <- try(gamlss(y~1, sigma.fo=~1, family="NEE",
                control=gamlss.control(n.cyc=2500, trace=FALSE)))
  
  if (class(mod)[1] == "try-error")
    res <- rep(NA, 2)
  else
    res <- c(exp(coef(mod, what="mu")), exp(coef(mod, what="sigma")))
  res
}

# Super function to simulate and write the estimated parameters
simul <- function(n) {
  result <- t(replicate(n=nrep, expr=simul_one(size=n)))
  result <- cbind(result, n)
  write(x=t(result), file="Simulations/simul_without_cov.txt", 
        ncol=3, append=TRUE)
}

# Code to generate the simulations given n --------------------------------

# Aqui se definen los valores de tamano muestral n
# Luego se define el numero de repeticiones
n <- seq(from=100, to=1000, by=100)
nrep <- 15

values <- expand.grid(n=n)
values
apply(values, 1, simul)


# Plots -------------------------------------------------------------------
dt <- read.table("Simulations/simul_without_cov.txt", 
                 col.names=c("mu_hat", "sigma_hat", "n"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# Numero de observaciones por cada n
num <- dt %>% group_by(n) %>% count()
mean(num$nn)
min(num$nn)

# Para obtener la metricas
trim <- 0.05

res <- dt %>% 
  drop_na() %>% 
  group_by(n) %>% 
  summarise(nobs=n(),
            mean_mu=mean(mu_hat, trim=trim), 
            mean_sigma=mean(sigma_hat, trim=trim),
            bias_mu=true_mu - mean(mu_hat, trim=trim), 
            bias_sigma=true_sigma - mean(sigma_hat, trim=trim),
            mse_mu=mean((true_mu - mu_hat)^2, trim=trim), 
            mse_sigma=mean((true_sigma - sigma_hat)^2), trim=trim)

res

# Mean -----------------------------------------------------
p1 <- ggplot(data=res, aes(x=n, y=mean_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(mu))) +
  geom_line(y=true_mu, col="red", lty="dashed")

p2 <- ggplot(data=res, aes(x=n, y=mean_sigma)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(sigma))) +
  geom_line(y=true_sigma, col="red", lty="dashed")

mean1 <- grid.arrange(p1, p2, nrow = 1)
mean1
ggsave(filename="Figs/mean1.pdf", 
       plot=mean1, 
       width=10, height=4)

# MSE -----------------------------------------------------
p1 <- ggplot(data=res, aes(x=n, y=mse_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(mu)))

p2 <- ggplot(data=res, aes(x=n, y=mse_sigma)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(sigma)))

mse1 <- grid.arrange(p1, p2, nrow = 1)
mse1
ggsave(filename="Figs/mse1.pdf", 
       plot=mse1, 
       width=10, height=4)

