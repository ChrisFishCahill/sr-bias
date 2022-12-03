#-------------------------------------------------------------------------------
# simulate su and peterson
# cahill, walters, punt, november 2022

# objectives:
# 1) simulate the su and peterson state space ricker model
# 2) estimate it in Stan--complete a self test to show model works
# 3) repeat for different depletion levels and exploitation trajectories

#-------------------------------------------------------------------------------
# packages
#-------------------------------------------------------------------------------
library(rstan)
library(tidybayes)
library(tidyverse)
library(cowplot)

# devtools::install_github("ChrisFishCahill/gg-qfc")
library(ggqfc)
library(future)
library(furrr)

#-------------------------------------------------------------------------------
# function to get su peterson dynamics
#-------------------------------------------------------------------------------

sr_model <- function(Ut = NA) {
  wt <- rnorm(n_year - k, 0, sdp) # process noise
  vt <- rnorm(n_year, 0, sdo) # observation noise

  # Initialize S, R, C
  S[1:k] <- R[1:k] <- ar / b # S' = R' = ln(a)/b = equilibrium S and R
  C[1:k] <- Ut[1:k] * R[1:k]

  # sequentially generate new recruits, catch, and spawners
  for (t in 1:(n_year - k)) {
    R[t + k] <- a * S[t] * exp(-b * S[t] + wt[t]) # truth + process noise
    C[t + k] <- Ut[t + k] * R[t + k]
    S[t + k] <- R[t + k] - C[t + k]
  }

  E <- S * exp(vt) # add obs. noise to S to get E

  out <- tibble(
    "E" = E,
    "R" = R,
    "S" = S,
    "ln_RS" = c(
      log(R[(k + 1):n_year] / S[1:(n_year - k)]),
      rep(NA, k)
    ),
    "C" = C,
    "Ut" = Ut,
    "wt" = c(rep(NA, k), wt),
    "vt" = vt,
    "year" = 1:n_year
  )
  out
}

#-------------------------------------------------------------------------------
# now simulate/estimate many times
#-------------------------------------------------------------------------------

get_fit <- function(sim = NA, Umax = NA, 
                    scenario = c("depleted", "recovering", "declining")
                    ) {
  # set up exploitation rate sequence
  if (scenario == "depleted") {
    # fish to low state determined by Umax
    Ut <- rep(NA, n_year)
    relU <- seq(from = 0, to = 1, by = 0.05)
    Ut[1:length(relU)] <- relU
    Ut[which(is.na(Ut))] <- 1
    Ut <- Ut * Umax
  } else if (scenario == "recovering") {
    # fish to low state, and then reduce exploitation 
    Ut <- rep(NA, n_year)
    relseq <- seq(from = 0, to = 1, by = 0.05)
    relU <- c(relseq, rep(1, n_year/2 - length(relseq))) * 0.585 
    # ramp Ut down to Umax: 
    relU <- c(relU, ifelse(rev(relseq)*0.585 > Umax, rev(relseq)*0.585, Umax))  
    Ut[1:length(relU)] <- relU
    Ut[which(is.na(Ut))] <- Umax
  } else if (scenario == "declining") {
    # almost no exploitation until year fifty, then increasing exploitation
    Ut <- rep(NA, n_year)
    Ut[1:(n_year / 2)] <- 0.01 
    relU <- seq(from = 0.01, to = 1, by = 0.05)
    relU <- c(relU, rep(1, (length(Ut) - sum(is.na(Ut)) - length(relU))))
    Ut[which(is.na(Ut))] <- c(relU[1], relU[2:length(relU)]*Umax) 
  }

  sim_dat <- sr_model(Ut = Ut) # draw a single time series|Ut

  # ---------------------------------------
  # take last n years to illustrate ts bias
  # ---------------------------------------
  n <- n_year - 50
  E <- sim_dat$E[(n_year - (n - 1)):n_year]
  C <- sim_dat$C[(n_year - (n - 3)):n_year]
  S <- sim_dat$S[(n_year - (n - 1)):n_year]

  inits <- function() {
    list(
      "ar" = ar,
      "ln_So" = rep(log(sim_dat$S[1], k)),
      "sdo" = sdo,
      "sdp" = sdp,
      "R" = sim_dat$R[(k + 1):n_year]
    )
  }
  if(scenario == "recovering") inits <- "random"
  stan_data <-
    list(
      "k" = k,
      "n_year" = length(E),
      "E" = E,
      "C" = C,
      ar_prior = c(ar, 0.25),
      ln_sdp_prior = c(log(sdp), 0.15),
      ln_sdo_prior = c(log(sdo), 0.15)
    )

  fit <-
    rstan::sampling(
      m,
      data = stan_data,
      #init = inits,
      pars = c("ar", "b", "sdo", "sdp", "smsy", "hmsy"),
      iter = 5000, warmup = 2500, chains = 1,
      control = list(adapt_delta = 0.999, max_treedepth = 15),
      verbose = FALSE
    )
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  divergent <- sapply(sampler_params, function(x) max(x[, "divergent__"]))
  res <- fit %>%
    spread_draws(ar, b, sdo, sdp, smsy, hmsy) %>%
    summarise_draws()
  res <- res %>% add_column(
    sim = sim, Umax = Umax, n_year = n,
    Dbar = mean(S) / sim_dat$S[1], 
    divergent = divergent, scenario = scenario
  )
  res
}

#-------------------------------------------------------------------------------
# set up some leading parameters / values for the f(x)
#-------------------------------------------------------------------------------

k <- 2 # age at maturity
n_year <- 100
ar <- b <- 1 # ln( ricker alpha), ricker b
hmsy <- 0.5 * ar - 0.07 * (ar^2) # su and peterson relationships
smsy <- hmsy / b # su and peterson relationships
a <- exp(ar)
sdp <- 0.001 # process error sd
sdo <- 0.15 # observation error sd-- note this affects catch, not calculation of R and S
E <- S <- rep(NA, n_year) # Escapement, Stock
C <- R <- rep(NA, n_year) # Catch, Recruits

# fish to low state determined by Umax
Umax <- 0.5999
Ut <- rep(NA, n_year)
relU <- seq(from = 0, to = 1, by = 0.05)
Ut[1:length(relU)] <- relU
Ut[which(is.na(Ut))] <- 1
Ut <- Ut * Umax

# visualize for jim:
sim_dat <- sr_model(Ut = Ut) 

sim_dat %>%
  ggplot(aes(x = S, y = R))+
  geom_point() + 
  theme_qfc()

# or not using ggplot:
plot(sim_dat$R ~ sim_dat$S, xlab = "stock", ylab = "recruits")

#-------------------------------------------------------------------------------
# call stan, fit the models with furrr
#-------------------------------------------------------------------------------
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

# compile the stan model
path <- "src/ss_ricker.stan"
m <- rstan::stan_model(path, verbose = T)

#  debugging
# set.seed(1)
# out = get_fit(sim = 1)
# out = pivot_longer(out, variable)
# system.time(
# out <- purrr::pmap_dfr(to_sim, get_fit)
# )

# exploitation rate maximums to simulate across
Umax <- seq(from = 0.01, to = .5999, length.out = 6)

scenario = c("depleted", "recovering", "declining")
sim <- seq_len(10)
to_sim <- expand.grid(sim = sim, Umax = Umax, scenario = scenario)

future::plan(multisession)

system.time({
  out <- future_pmap_dfr(to_sim, get_fit,
    .options = furrr_options(seed = TRUE, 
                             scheduling = Inf),
    .progress = TRUE
  )
})

#saveRDS(out, file = "sims/su-peterson-sims.rds")
summary(warnings())

#-------------------------------------------------------------------------------
# testing--call simulation and estimate model one time
#-------------------------------------------------------------------------------

# set.seed(3)
# dat <- sr_model()
#
# options(mc.cores = parallel::detectCores())
# rstan::rstan_options(auto_write = TRUE)
#
# # compile the stan model
# path <- "src/ss_ricker.stan"
# m <- rstan::stan_model(path, verbose = T)
#

# Different exploitation scenarios

# depleted low state scenario:
# Ut <- rep(NA, n_year)
# Umax <- .59
# relU <- seq(from = 0, to = 1, by = 0.05)
# Ut[1:length(relU)] <- relU
# Ut[which(is.na(Ut))] <- 1
# Ut <- Ut * Umax
#
# dat <- sr_model(Ut = Ut)
# plot(dat$S, type = "b", ylab="Stock", xlab = "Year")
# plot(Ut)
#
# # recovering from depleted state
# Ut <- rep(NA, n_year)
# Umax <- .59
# relseq <- seq(from = 0, to = 1, by = 0.05)
# relU <- c(relseq, rep(1, n_year/2 - length(relseq))) * 0.59 # first half
# #
# # # ramp down to Umax:
# relU <- c(relU, ifelse(rev(relseq)*0.59 > Umax, rev(relseq)*0.59, Umax))
# Ut[1:length(relU)] <- relU
# Ut[which(is.na(Ut))] <- Umax
#
# dat <- sr_model(Ut = Ut)
# par(mfrow=c(2,1))
# plot(Ut, type= "b", ylab="Ut", xlab="Year")
# plot(dat$S, type = "b", ylab="Stock", xlab = "Year")
#
# # declining from pristine state
# Umax <- 0.59
# Ut <- rep(NA, n_year)
# Ut[1:(n_year / 2)] <- 0.01 # fish heavily first half of sim
# relU <- seq(from = 0.01, to = 1, by = 0.05)
# relU <- c(relU, rep(1, (length(Ut) - sum(is.na(Ut)) - length(relU))))
# Ut[which(is.na(Ut))] <- c(relU[1], relU[2:length(relU)]*Umax) # fish at Umax
#
# plot(Ut, type = "b")
# dat <- sr_model(Ut = Ut)
# plot(dat$S, type = "b", ylab="Stock", xlab = "Year")
#
# # set.seed(3)
# # recovering from depleted state
# Ut <- rep(NA, n_year)
# Umax <- .4
# relseq <- seq(from = 0, to = 1, by = 0.05)
# relU <- c(relseq, rep(1, n_year/2 - length(relseq))) * 0.59 # first half
# #
# # # ramp down to Umax:
# relU <- c(relU, ifelse(rev(relseq)*0.59 > Umax, rev(relseq)*0.59, Umax))
# Ut[1:length(relU)] <- relU
# Ut[which(is.na(Ut))] <- Umax
# dat <- sr_model(Ut = Ut)
# plot(dat$S, type = "b")
#
# # take last n years to illustrate the time series bias
# n <- 50
# E <- dat$E[(n_year - (n - 1)):n_year]
# C <- dat$C[(n_year - (n - 3)):n_year]
#
# inits <- function() {
#   list(
#     "ar" = ar,
#     "ln_So" = rep(log(dat$S[1], k)),
#     "sdo" = sdo,
#     "sdp" = sdp,
#     "R" = dat$R[(k + 1):n_year]
#   )
# }
#
# stan_data <-
#   list(
#     "k" = k,
#     "n_year" = length(E),
#     "E" = E,
#     "C" = C,
#     ar_prior = c(ar, 0.25),
#     ln_sdp_prior = c(log(sdp), 0.15),
#     ln_sdo_prior = c(log(sdo), 0.15)
#   )
#
# fit <-
#   rstan::sampling(
#     m,
#     data = stan_data,
#     pars = c("ar", "b", "sdo", "sdp", "smsy", "hmsy", "R", "S"),
#     init = inits,
#     iter = 100, warmup = 0, chains = 1,
#     control = list(adapt_delta = 0.999, max_treedepth = 13)
#   )

# pairs(fit, pars = c("ar", "b", "sdo", "sdp"))
# fit

# res = fit %>% spread_draws(ar, b, sdo, sdp) %>%
#   summarise_draws()
# shinystan::launch_shinystan(fit)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# End end end