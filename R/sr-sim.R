#-------------------------------------------------------------------------------
# simulate su and peterson
# cahill, punt, november 2022

# objectives:
# 1) simulate the su and peterson state space ricker model
# 2) estimate it in Stan--complete a self test to show the model works
# 3) Reduce data quality from the simulated trajectory of catches
#    to show time-series bias
# 4) visualize
#-------------------------------------------------------------------------------
# packages 
#-------------------------------------------------------------------------------
library(tidybayes)
library(tidyverse)
library(cowplot)
library(rstan)
library(tidybayes)
library(cowplot)
# devtools::install_github("ChrisFishCahill/gg-qfc")
library(ggqfc)

#-------------------------------------------------------------------------------
# function to get su peterson dynamics
#-------------------------------------------------------------------------------

sr_model <- function() {
  wt <- rnorm(n_year - k, 0, sdp) # process noise
  vt <- rnorm(n_year, 0, sdo) # observation noise

  # Initialize S, R, C
  S[1:k] <- R[1:k] <- ar / b # S' = R' = ln(a)/b = equilibrium S
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
  return(out)
}

#-------------------------------------------------------------------------------
# set up some leading parameters / values for the f(x)
#-------------------------------------------------------------------------------

k <- 2 # age at maturity
n_year <- 100 # 50
ar <- b <- 1 # ln( ricker alpha), ricker b
a <- exp(ar)
sdp <- 0.05 # process error sd
sdo <- 0.3 # observation error sd
E <- S <- rep(NA, n_year) # Escapement, Stock
C <- R <- rep(NA, n_year) # Catch, Recruits

# set up exploitation rate sequence
Ut <- rep(NA, n_year)
U <- 0.6
relU <- seq(from = 0, to = 1, by = 0.05)
Ut[1:length(relU)] <- relU
Ut[which(is.na(Ut))] <- 1
Ut <- Ut * U

#--------------------------------------------------------------------------
# call the f(x) and estimate it once in stan
#--------------------------------------------------------------------------
set.seed(3)
dat <- sr_model()

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

# compile the stan model
path <- "src/ss_ricker.stan"
m <- rstan::stan_model(path, verbose = T)

set.seed(3)
dat <- sr_model()

# take last n years to illustrate the time series bias
n <- 50
E <- dat$E[(n_year - (n - 1)):n_year]
C <- dat$C[(n_year - (n - 3)):n_year]

inits <- function() {
  list(
    "ar" = ar,
    "ln_So" = rep(log(dat$S[1], k)),
    "sdo" = sdo,
    "sdp" = sdp,
    "R" = dat$R[(k + 1):n_year]
  )
}

stan_data <-
  list(
    "k" = k,
    "n_year" = length(E),
    "E" = E,
    "C" = C,
    ar_prior = c(ar, 0.2),
    ln_sdp_prior = c(log(sdp), 0.15),
    ln_sdo_prior = c(log(sdo), 0.15)
  )

fit <-
  rstan::sampling(
    m,
    data = stan_data,
    init = inits,
    pars = c("ar", "b", "sdo", "sdp", "R", "S"),
    iter = 5000, warmup = 2500, chains = 1
  )

fit
# res = fit %>% spread_draws(ar, b, sdo, sdp) %>%
#   summarise_draws()

#-------------------------------------------------------------------------------
# now simulate/estimate many times
#-------------------------------------------------------------------------------

get_fit <- function(sim = NA) {
  sim_dat <- sr_model() # draw a single time series from 1:n_year

  # ---------------------------------------
  # fit to all years from simulation
  # ---------------------------------------
  n <- n_year
  E <- sim_dat$E[(n_year - (n - 1)):n_year]
  C <- sim_dat$C[(n_year - (n - 3)):n_year]

  inits <- function() {
    list(
      "ar" = ar,
      "ln_So" = rep(log(sim_dat$S[1], k)),
      "sdo" = sdo,
      "sdp" = sdp,
      "R" = sim_dat$R[(k + 1):n_year]
    )
  }

  stan_data <-
    list(
      "k" = k,
      "n_year" = length(E),
      "E" = E,
      "C" = C,
      ar_prior = c(ar, 0.2),
      ln_sdp_prior = c(log(sdp), 0.15),
      ln_sdo_prior = c(log(sdo), 0.15)
    )

  fit <-
    rstan::sampling(
      m,
      data = stan_data,
      init = inits,
      pars = c("ar", "b", "sdo", "sdp"),
      iter = 5000, warmup = 2500, chains = 1,
      control = list(adapt_delta = 0.99, max_treedepth = 13)
    )

  res <- fit %>%
    spread_draws(ar, b, sdo, sdp) %>%
    summarise_draws()

  res <- res %>% add_column(sim = sim, n_year = n_year)
  # ---------------------------------------
  # take last n years to illustrate ts bias
  # ---------------------------------------
  n <- 50
  E <- sim_dat$E[(n_year - (n - 1)):n_year]
  C <- sim_dat$C[(n_year - (n - 3)):n_year]

  inits <- function() {
    list(
      "ar" = ar,
      "ln_So" = rep(log(sim_dat$S[1], k)),
      "sdo" = sdo,
      "sdp" = sdp,
      "R" = sim_dat$R[(k + 1):n_year]
    )
  }

  stan_data <-
    list(
      "k" = k,
      "n_year" = length(E),
      "E" = E,
      "C" = C,
      ar_prior = c(ar, 0.2),
      ln_sdp_prior = c(log(sdp), 0.15),
      ln_sdo_prior = c(log(sdo), 0.15)
    )

  fit <-
    rstan::sampling(
      m,
      data = stan_data,
      init = inits,
      pars = c("ar", "b", "sdo", "sdp"),
      iter = 5000, warmup = 2500, chains = 1,
      control = list(adapt_delta = 0.99, max_treedepth = 13)
    )

  res2 <- fit %>%
    spread_draws(ar, b, sdo, sdp) %>%
    summarise_draws()
  res2 <- res2 %>% add_column(sim = sim, n_year = 50)
  res <- bind_rows(res, res2)
  res
}

#  testing/debugging
# set.seed(1)
# out = get_fit(sim = 1)
# out = pivot_longer(out, variable)
# system.time(
#   out <- purrr::pmap_dfr(to_sim, get_fit)
# )

to_sim <- tibble(sim = seq_len(100))

future::plan(multisession)

system.time({
  out <- future_pmap_dfr(to_sim, get_fit,
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )
})
out <- out %>% pivot_longer(variable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# plot the output from simulation

ap <- out %>%
  filter(value == "ar") %>%
  ggplot(aes(y = median, x = as.factor(n_year))) +
  geom_violin(width = 0.75) +
  geom_jitter(width = 0.12, alpha = 0.5) +
  geom_hline(yintercept = c(ar), lty = 2, color = "steelblue", lwd = 2) +
  ggtitle(expression(ln ~ alpha)) +
  ylab("Posterior median estimates") +
  xlab("Low or high data quality (final 50 or all 100 yrs of sim)") +
  theme_qfc() +
  stat_summary(fun = median, geom = "point", size = 3, col = "darkorange3")
ap

bp <- out %>%
  filter(value == "b") %>%
  ggplot(aes(y = median, x = as.factor(n_year))) +
  geom_violin(width = 1.5) +
  geom_jitter(width = 0.05, alpha = 0.5) +
  geom_hline(yintercept = b, lty = 2, color = "steelblue", lwd = 2) +
  ggtitle(expression(beta)) +
  ylab("Posterior median estimates") +
  xlab("Low or high data quality (final 50 or all 100 yrs of sim)") +
  theme_qfc() +
  stat_summary(fun = median, geom = "point", size = 3, col = "darkorange3")
bp

sigo <- out %>%
  filter(value == "sdo") %>%
  ggplot(aes(y = median, x = as.factor(n_year))) +
  geom_violin(width = 0.12) +
  geom_jitter(width = 0.05, alpha = 0.5) +
  geom_hline(yintercept = sdo, lty = 2, color = "steelblue", lwd = 2) +
  ggtitle(expression(sigma[observation])) +
  ylab("Posterior median estimates") +
  xlab("Low or high data quality (final 50 or all 100 yrs of sim)") +
  theme_qfc() +
  stat_summary(fun = median, geom = "point", size = 3, col = "darkorange3")
sigo

sigp <- out %>%
  filter(value == "sdp") %>%
  ggplot(aes(y = median, x = as.factor(n_year))) +
  geom_violin(width = 0.12) +
  geom_jitter(width = 0.05, alpha = 0.5) +
  geom_hline(yintercept = sdp, lty = 2, color = "steelblue", lwd = 2) +
  ggtitle(expression(sigma[process])) +
  ylab("Posterior median estimates") +
  xlab("Low or high data quality (final 50 or all 100 yrs of sim)") +
  theme_qfc() +
  stat_summary(fun = median, geom = "point", size = 3, col = "darkorange3")
sigp

p <- plot_grid(ap, bp, sigo, sigp, ncol = 2)
p

ggsave("plots/ts-bias-demonstration.pdf", width = 11, height = 8)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# call the function and make some more plots
set.seed(3)
dat <- sr_model()

p1 <- dat %>%
  ggplot(aes(x = year, y = S)) +
  geom_point() +
  geom_line() +
  theme_qfc() +
  ylim(0, 1.0) +
  ggtitle("True S vs. time") +
  geom_hline(yintercept = 0.1, linetype = 2)

p2 <- dat %>%
  ggplot(aes(x = year, y = Ut)) +
  geom_point() +
  geom_line() +
  theme_qfc() +
  ylab("Exploitation rate U(t)") +
  ggtitle("Ut vs. time")

p3 <- dat %>%
  ggplot(aes(x = S, y = ln_RS)) +
  geom_point() +
  ylab("ln(R/S)") +
  xlab("S") +
  ggtitle("True ln(R/S) vs. S relationship") +
  theme_qfc()

p4 <- dat %>%
  ggplot(aes(x = year, y = E)) +
  geom_point() +
  xlab("year") +
  ylab("Escapement") +
  theme_qfc() +
  geom_hline(yintercept = 0.1, linetype = 2) +
  ggtitle("Observed Escapement vs. time")

p5 <- dat %>%
  ggplot(aes(x = year, y = wt)) +
  geom_line() +
  xlab("year") +
  ylab(expression(sigma[process])) +
  theme_qfc() +
  ggtitle(expression(sigma[process] ~ vs. ~ time))

p6 <- dat %>%
  ggplot(aes(x = year, y = vt)) +
  geom_line() +
  xlab("year") +
  ylab(expression(sigma[obs])) +
  theme_qfc() +
  ggtitle(expression(sigma[obs] ~ vs. ~ time))

p <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)
p

ggsave("plots/ts-simulation-demonstration.pdf", width = 11, height = 8)


# next chunk requires a single fit called fit:

R2 <- fit %>% spread_draws(R[year])
R2 <- R2 %>%
  filter(year > k) %>%
  mutate(year = year - k)
S2 <- fit %>% spread_draws(S[year])
S2 <- S2 %>% filter(year <= (n_year - k))

dat <- left_join(R2, S2)
dat$ln_RS <- log(dat$R / dat$S)

# pluck out the best estimated ar, b
ests <- fit %>%
  gather_draws(ar, b) %>%
  median_qi()

ar_est <- ests$.value[which(ests$.variable == "ar")]
b_est <- ests$.value[which(ests$.variable == "b")]

devs <- fit %>%
  spread_draws(ar, b) %>%
  sample_n(size = 300)

p4 <- dat %>%
  median_qi() %>%
  ggplot(aes(x = S, y = ln_RS, ymin = ln_RS.lower, ymax = ln_RS.upper)) +
  geom_pointinterval(alpha = 0.25) +
  ylab("ln(R/S)") +
  ggqfc::theme_qfc()

p4 <- p4 + geom_abline(
  intercept = devs$ar, slope = -devs$b, color = "black",
  linetype = 1, size = 0.5, alpha = .05
)

p4 <- p4 + geom_abline(
  intercept = ar_est, slope = -b_est, color = "black",
  linetype = 2, size = 1.5
)

p4 <- p4 + geom_abline(
  intercept = ar, slope = -b, color = "steelblue",
  linetype = 2, size = 2
)

p4 <- p4 + ggtitle("Estimated (black/gray) ln(R/S) vs. S  vs. true relationship (blue)")

set.seed(3)
dat2 <- sr_model()
dat2$color <- c(rep("not included", n_year / 2), rep("included", n_year / 2))
p5 <- dat2 %>%
  ggplot(aes(x = S, y = ln_RS, color = color)) +
  geom_point() +
  scale_color_manual(values = c("black", "steelblue")) +
  ylab("ln(R/S)") +
  xlab("S") +
  ggtitle("True ln(R/S) vs. S relationship") +
  theme_qfc()

p5

p6 <- dat2 %>%
  ggplot(aes(y = S, x = year, color = color)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("black", "steelblue")) +
  ylab("Stock Size (S)") +
  xlab("Year") +
  ggtitle("Visualizing what the assessment can and cannot 'see' in time series bias scenario") +
  theme_qfc()
p6

p <- plot_grid(p6, p5, p4, ncol = 1)

ggsave("plots/ts-bias-demonstration-one-fit.pdf", width = 8, height = 11)

#-------------------------------------------------------------------------------
# even more plotting
#-------------------------------------------------------------------------------
#
# p1 <- fit %>%
#   spread_draws(R[year]) %>%
#   median_qi() %>%
#   ggplot(aes(x = year, y = R, ymin = .lower, ymax = .upper)) +
#   geom_pointinterval()
# p1
# dat <- data.frame(year = 1:(n - k), Rtrue = R[(n_year - (n - (k + 1))):n_year])
# p1 <- p1 + geom_point(
#   data = dat, aes(
#     x = year, y = Rtrue, ymin = Rtrue,
#     ymax = Rtrue
#   ), shape = 16, color = "red",
#   size = 2
# ) +
#   ggqfc::theme_qfc()
#
# p1
#
# p2 <- fit %>%
#   spread_draws(S[year]) %>%
#   median_qi() %>%
#   ggplot(aes(x = year, y = S, ymin = .lower, ymax = .upper)) +
#   geom_pointinterval()
# p2
# dat <- data.frame(year = 1:n, Strue = S[(n_year - (n - 1)):n_year])
# p2 <- p2 + geom_point(
#   data = dat, aes(
#     x = year, y = Strue, ymin = Strue,
#     ymax = Strue
#   ), shape = 16, color = "red",
#   size = 2
# ) +
#   ggqfc::theme_qfc()
# p2
#
# p3 <- fit %>%
#   gather_draws(ar, b, sdp, sdo) %>%
#   median_qi() %>%
#   ggplot(aes(x = .variable, ymin = .lower, ymax = .upper, y = .value)) +
#   geom_pointinterval() +
#   xlab("Parameter") +
#   ylab("Value") +
#   ggqfc::theme_qfc()
# dat <- data.frame(
#   .variable = c("ar", "b", "sdo", "sdp"),
#   .value = c(ar, b, sdo, sdp)
# )
# p3 <- p3 + geom_point(
#   data = dat, aes(
#     x = .variable, y = .value, ymin = .value,
#     ymax = .value
#   ),
#   shape = 16, color = "red",
#   size = 2
# )
# p3
#
# p <- plot_grid(p3, p1, p2, ncol = 1)
# p
# ggsave("plots/self-test-su-peterson.pdf", width = 8, height = 11)

# create a few posteriors for stock-recruit, i.e., ln(R/S) vs. S