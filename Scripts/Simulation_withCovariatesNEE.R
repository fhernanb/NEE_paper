library(gamlss)
library(RelDists)

#-----------------------------------------------------------------
#------------------ EJEMPLO CON LA DISTRIBUCION ------------------
#----------------------- Zero inflated Poisson -------------------
#---------------------------- NEE --------------------------------
#-----------------------------------------------------------------

# En la distribucion NEE
# mu > 0,         por tanto usaremos funcion de enlace log
# sigma > 0  por tanto usaremos funcion de enlace log

# Creando las funciones de enlace inversas
#logit_inv <- function(x) exp(x) / (1 + exp(x))
# No vamos a crear log_inv porque esa funcion ya existe y 
# llama exp( )

# The parameters ----------------------------------------------------------
# Los siguientes son los valores de los betas para el modelo de regresion
# debemos chequear que esos numeros den valores correctos de mu y sigma
# no podemos asignar numeros a la loca

true_b0_mu <- -1.4    # intercept for mu
true_b1_mu <-  4.6   # slope for mu
true_g0_si <-  2.1    # intercept for sigma
true_g1_si <-  2.3   # slope for sigma

# Useful functions to the simulation study --------------------------------

# Funcion para obtener mu_hat y sigma_hat para un valor fijo de n
simul_one <- function(size) {
  x1 <- runif(n=size)
  x2 <- runif(n=size)
  mu    <- exp(true_b0_mu + true_b1_mu * x1) # Aprox mu = 2
  sigma <- exp(true_g0_si + true_g1_si * x2) # Aprox sigma = 255
  y <- rNEE(n=size, mu=mu, sigma=sigma)
  mod <- NULL
  mod <- try(gamlss(y~x1, sigma.fo=~x2, family="NEE",
                    control=gamlss.control(n.cyc=2500, trace=FALSE)))
  if (class(mod)[1] == "try-error")
    res <- rep(NA, 4)
  else
    res <- c(coef(mod, what="mu"), coef(mod, what="sigma"))
  res
}

# Super function to simulate and write the estimated parameters
simul <- function(n) {
  result <- t(replicate(n=nrep, expr=simul_one(size=n)))
  result <- cbind(result, n)
  write(x=t(result), file="Simulations/simul_with_cov.txt", 
        ncol=5, append=TRUE)
}

# Code to generate the simulations given n --------------------------------

# Aqui se definen los valores de tamano muestral n
# Luego se define el numero de repeticiones
n <- seq(from=100, to=1000, by=100)
nrep <- 90

values <- expand.grid(n=n)
values
apply(values, 1, simul)

# Plots -------------------------------------------------------------------
dt <- read.table("Simulations/simul_with_cov.txt", 
                 col.names=c("b0_mu_hat", "b1_mu_hat", 
                             "g0_si_hat", "g1_si_hat", 
                             "n"))

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
            mean_b0_mu=mean(b0_mu_hat, trim=trim),
            mean_b1_mu=mean(b1_mu_hat, trim=trim),
            mean_b0_si=mean(g0_si_hat, trim=trim),
            mean_b1_si=mean(g1_si_hat, trim=trim),
            mse_b0_mu=mean((true_b0_mu - b0_mu_hat)^2, trim=trim), 
            mse_b1_mu=mean((true_b1_mu - b1_mu_hat)^2, trim=trim),
            mse_b0_si=mean((true_g0_si - g0_si_hat)^2, trim=trim), 
            mse_b1_si=mean((true_g1_si - g1_si_hat)^2, trim=trim))

res

# Mean -----------------------------------------------------
p1 <- ggplot(data=res, aes(x=n, y=mean_b0_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(beta)[0]), 
       title=expression("Mean for the intercept in" ~ mu)) +
  geom_line(y=true_b0_mu, col="red", lty="dashed")

p2 <- ggplot(data=res, aes(x=n, y=mean_b1_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(beta)[1]),
       title=expression("Mean for the slope in" ~ mu)) +
  geom_line(y=true_b1_mu, col="red", lty="dashed")

p3 <- ggplot(data=res, aes(x=n, y=mean_b0_si)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(gamma)[0]),
       title=expression("Mean for the intercept in" ~ sigma)) +
  geom_line(y=true_g0_si, col="red", lty="dashed")

p4 <- ggplot(data=res, aes(x=n, y=mean_b1_si)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(gamma)[1]),
       title=expression("Mean for the slope in" ~ sigma)) +
  geom_line(y=true_g1_si, col="red", lty="dashed")

mean2 <- grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
mean2
ggsave(filename="mean2.pdf", 
       plot=mean2, 
       width=10, height=8)

# MSE -----------------------------------------------------
p1 <- ggplot(data=res, aes(x=n, y=mse_b0_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(beta)[0]),
       title=expression("MSE for the intercept in" ~ mu))

p2 <- ggplot(data=res, aes(x=n, y=mse_b1_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(beta)[1]),
       title=expression("MSE for the slope in" ~ mu))

p3 <- ggplot(data=res, aes(x=n, y=mse_b0_si)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(gamma)[0]),
       title=expression("MSE for the intercept in" ~ sigma))

p4 <- ggplot(data=res, aes(x=n, y=mse_b1_si)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(gamma)[1]),
       title=expression("MSE for the slope in" ~ sigma))

mse2 <- grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
mse2
ggsave(filename="mse2.pdf", 
       plot=mse2, 
       width=10, height=8)

