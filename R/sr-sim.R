# simulate su and peterson
# cahill, punt november 2022

k <- 2 # age at maturity
n_year <- 50
ar <- 1 # ln( ricker alpha)
a <- exp(ar)
b <- 1
sdp <- 0.05 # process error sd
sdo <- 0.3 # observation error sd

E <- S <- rep(NA, n_year) # Escapement, Stock
C <- R <- rep(NA, n_year) # Catch, Recruits

# set up exploitation rate sequence
Ut <- rep(NA, n_year)
U <- 0.5
relU <- seq(from = 0, to = 1, by = 0.02)
Ut[1:length(relU)] <- relU
Ut[which(is.na(Ut))] <- 1
Ut <- Ut * U

set.seed(1)
wt <- rnorm(n_year - k, 0, sd = sdp) # process noise
# Initialize S, R, C
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
library(tidyverse)
library(tidybayes)

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
    "C" = C[(k + 1):n_year],
    ar_prior = c(1, 3),
    b_prior = c(1, 3),
    sdp_prior = c(0, 0.5),
    sdo_prior = c(0, 1),
    So_prior = c(0, 3)
  )

inits <- function() {
  list(
    "ar" = jitter(ar),
    "b" = jitter(b),
    "sdo" = jitter(log(sdo)),
    "sdp" = jitter(log(sdp)),
    "R" = rep(jitter(1), length((k + 1):n_year)),
    "So" = jitter(c(1, 1))
  )
}

fit <-
  rstan::sampling(
    m,
    data = stan_data,
    pars = c("ar", "b", "sdo", "sdp", "So", "R", "S"),
    iter = 3000, chains = 4,
    control = list(adapt_delta = 0.999, max_treedepth = 12)
  )

# shinystan::launch_shinystan(fit)

#-------------------------------------------------------------------------------
# plot stuff
library(tidybayes)
library(tidyverse)
library(cowplot)

# see how we stack up against truth
ar
b
sdo
sdp
S[1:2] # So 1,2

# examine the fit
fit

p1 <- fit %>%
  spread_draws(R[year]) %>%
  median_qi() %>%
  ggplot(aes(x = year, y = R, ymin = .lower, ymax = .upper)) +
  geom_pointinterval()
p1
dat <- data.frame(year = 1:(n_year - k), Rtrue = R[(k + 1):n_year])
p1 <- p1 + geom_point(
  data = dat, aes(
    x = year, y = Rtrue, ymin = Rtrue,
    ymax = Rtrue
  ), shape = 16, color = "red",
  size = 2
) +
  ggqfc::theme_qfc()

p1

p2 <- fit %>%
  spread_draws(S[year]) %>%
  median_qi() %>%
  ggplot(aes(x = year, y = S, ymin = .lower, ymax = .upper)) +
  geom_pointinterval()
p2
dat <- data.frame(year = 1:n_year, Strue = S)
p2 <- p2 + geom_point(
  data = dat, aes(
    x = year, y = Strue, ymin = Strue,
    ymax = Strue
  ), shape = 16, color = "red",
  size = 2
) +
  ggqfc::theme_qfc()
p2


p3 <- fit %>%
  gather_draws(ar, b, sdp, sdo) %>%
  median_qi() %>%
  ggplot(aes(x = .variable, ymin = .lower, ymax = .upper, y = .value)) +
  geom_pointinterval() +
  xlab("Parameter") +
  ylab("Value") +
  ggqfc::theme_qfc()
dat <- data.frame(
  .variable = c("ar", "b", "sdo", "sdp"),
  .value = c(ar, b, sdo, sdp)
)
p3 <- p3 + geom_point(
  data = dat, aes(
    x = .variable, y = .value, ymin = .value,
    ymax = .value
  ),
  shape = 16, color = "red",
  size = 2
)

p <- plot_grid(p3, p1, p2, ncol = 1)
p
ggsave("plots/self-test-su-peterson.pdf", width = 8, height = 11)

#-------------------------------------------------------------------------------
# now try TMB -- not yet working
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
  "wt" = rep(0, length((k + 1):n_year))
)

cppfile <- "src/ss_ricker.cpp"
compile(cppfile)
dyn.load(dynlib("src/ss_ricker"))
random <- "wt"
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
abline(0, 1)
plot(obj$report(opt$par)$`S` ~ S[1:n_year])
abline(0, 1)
obj$report()
