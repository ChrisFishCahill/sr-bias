# simulate su and peterson
# cahill, punt november 2022

k <- 2 # age at maturity
n_year <- 100
ar <- 1 # ln( ricker alpha)
a <- exp(ar) 
b <- 1
sdp <- 0.05 # process error sd
sdo <- 0.3 # observation error sd

E <- S <- rep(NA, n_year) # Escapement, Stock
C <- R <- h <- rep(NA, n_year) # Catch, Recruits, finite harvest rate

# Initialize S, R, C
Ut <- rep(NA, n_year)
U = 0.5
relU <- seq(from = 0, to = 1, by = 0.02)
Ut[1:length(relU)] <- relU
Ut[which(is.na(Ut))] <- 1
Ut <- Ut * U

#set.seed(1)
# h[(k + 1):n_year] <- runif(n_year - k, 0.1, 0.35) # harvest rate
wt <- rnorm(n_year-k, 0, sd = sdp) # process noise

S[1:k] <- ar / b # S' = ln(a)/b = equilibrium S
R[1:k] <- ar / b # R' = ln(a)/b = equilibrium R
C[1:k] <- Ut[1:k] * R[1:k] # equilibrium C given constant h = 0.1

# sequentially generate new recruits, catch, and spawners
for (t in 1:(n_year - k)) {
  R[t + k] <- a * S[t] * exp(-b * S[t] + wt[t]) # truth + process noise
  C[t + k] <- Ut[t + k] * R[t + k]
  S[t + k] <- R[t + k] - C[t + k]
}

# now add in observation noise
vt <- rnorm(n_year, 0, sdo)
E <- S * exp(vt) # obs. noise

# plot it
plot(log(R / S) ~ S,
     ylab = "ln(R/S) vs. S",
     xlab = "S"
)
plot(R ~ S, "ylab" = "Recruits", xlab = "Stock")
plot(R ~ E, "ylab" = "Recruits", xlab = "Escapement")

confint(lm(log(R / S) ~ S))
ar
b

#--------------------------------------------------------------------------
# estimate it in stan
#--------------------------------------------------------------------------
library(rstan)
options(mc.cores = parallel::detectCores())

# eliminate redundant compilations
rstan::rstan_options(auto_write = TRUE)

# compile the stan model
path <- "src/ss_ricker.stan"
m <- rstan::stan_model(path, verbose = T)

# set up the data and the initial values
stan_data <-
  list(
    "k" = k,
    "n_year" = n_year,
    "E" = E,
    "C" = C[(k+1):n_year], 
    ar_prior = c(1, 5), 
    b_prior = c(1, 5),
    sdp_prior = c(0, 1), 
    sdo_prior = c(0, 5), 
    So_prior = c(0, 5)
  )

inits <- function() {
  list(
    "ar" = jitter(ar),
    "b" = jitter(b),
    "sdo" = jitter(log(sdo)),
    "sdp" = jitter(log(sdp)),
    "R" = rep(jitter(1), length((k+1):n_year)), 
    "So" = jitter(c(1,1))
  )
}

fit <-
  rstan::sampling(
    m,
    data = stan_data,
    pars = c("ar", "b", "sdo", "sdp", "So", "R"),
    iter = 8000, chains = 4, 
    control=list(adapt_delta = 0.999, max_treedepth = 12)
  )

shinystan::launch_shinystan(fit)

library(tidybayes)
library(tidyverse)

fit %>%
  spread_draws(R[year]) %>%
  ggplot(aes(x = R)) + 
  geom_violin()

#-------------------------------------------------------------------------------
# now try TMB
#-------------------------------------------------------------------------------
library(TMB)

data <- list(
  "k" = k,
  "E" = E,
  "C" = C
)
par <- list(
  "ar" = ar,
  "br" = log(b),
  "ln_sdo" = log(sdo),
  "ln_sdp" = log(sdp),
  "wt" = rep(0, length((k+1):n_year))
)

cppfile <- "src/ss_ricker.cpp"
compile(cppfile)
dyn.load(dynlib("src/ss_ricker"))
random = "wt"
obj <- MakeADFun(data = data, parameters = par, random = random)

obj$fn(obj$par)
obj$gr(obj$par)
obj$report()

opt <- nlminb(
  start = obj$par, objective = obj$fn,
  gradient = obj$gr
)

opt <- nlminb(
  start = opt$par, objective = obj$fn,
  gradient = obj$gr
)
opt$convergence
opt$objective
obj$report()

opt$SD <- sdreport(obj)
opt$SD

plot(obj$report(opt$par)$`R` ~ R)
abline(0,1)
plot(obj$report(opt$par)$`S` ~ S[1:n_year])
abline(0,1)
obj$report()
