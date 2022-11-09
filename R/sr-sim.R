# simulate su and peterson
# cahill, punt, walters november 2022

library(tidybayes)
library(tidyverse)
library(cowplot)
library(rstan)
library(tidybayes)
library(cowplot)

# devtools::install_github("ChrisFishCahill/gg-qfc")
library(ggqfc)

k <- 2 # age at maturity
n_year <- 200 # 50
ar <- 1 # ln( ricker alpha)
a <- exp(ar)
b <- 1
sdp <- 0.05 # process error sd
sdo <- 0.3 # observation error sd

E <- S <- rep(NA, n_year) # Escapement, Stock
C <- R <- rep(NA, n_year) # Catch, Recruits

# set up exploitation rate sequence
Ut <- rep(NA, n_year)
U <- 0.6
relU <- seq(from = 0, to = 1, by = 0.025)
Ut[1:length(relU)] <- relU
Ut[which(is.na(Ut))] <- 1
Ut <- Ut * U

set.seed(56)
wt <- rnorm(n_year - k, 0, sd = sdp) # process noise
# Initialize S, R, C
S[1:k] <- ar / b # S' = ln(a)/b = equilibrium S
R[1:k] <- ar / b # R' = ln(a)/b = equilibrium R
C[1:k] <- Ut[1:k] * R[1:k]

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
plot(log(R[(k + 1):n_year] / S[1:(n_year - k)]) ~ S[1:(n_year - k)],
  ylab = "ln(R/S) vs. S",
  xlab = "S"
)
# plot(R ~ S, "ylab" = "Recruits", xlab = "Stock")
plot(R ~ E, "ylab" = "Recruits", xlab = "Escapement")

# confint(lm(log(R / S) ~ S))
# ar
# b

plot(S, xlab = "Year", type = "b")
abline(h = 0.1, lty = 3)

#--------------------------------------------------------------------------
# estimate it in stan
#--------------------------------------------------------------------------
options(mc.cores = parallel::detectCores())

# eliminate redundant compilations
rstan::rstan_options(auto_write = TRUE)

# compile the stan model
path <- "src/ss_ricker.stan"
m <- rstan::stan_model(path, verbose = T)

# set up the data and the initial values
inits <- function() {
  list(
    "ar" = ar,
    "ln_So" = rep(log(S[1], k)),
    "sdo" = jitter(sdo),
    "sdp" = jitter(sdp),
    "R" = R[(k + 1):n_year]
  )
}

# take last n years
n <- n_year
E2 <- E[(n_year - (n - 1)):n_year]
C2 <- C[(n_year - (n - 3)):n_year]

stan_data <-
  list(
    "k" = k,
    "n_year" = length(E2),
    "E" = E2,
    "C" = C2,
    # "E" = E, # E2
    # "C" = C[(k + 1):n_year],
    ar_prior = c(ar, 0.25),
    ln_So_prior = c(log(S[1]), 0.25),
    ln_sdp_prior = c(log(sdp), 0.15),
    ln_sdo_prior = c(log(sdo), 0.15)
  )

stan_data <-
  list(
    "k" = k,
    "n_year" = length(E2),
    "E" = E2,
    "C" = C2,
    # "E" = E, # E2
    # "C" = C[(k + 1):n_year],
    ar_prior = c(ar, 0.5),
    ln_So_prior = c(log(S[1]), 1),
    ln_sdp_prior = c(log(sdp), 0.5),
    ln_sdo_prior = c(log(sdo), 0.5)
  )

fit <-
  rstan::sampling(
    m,
    data = stan_data,
    pars = c("ar", "b", "sdo", "sdp", "R", "S"),
    iter = 5000, warmup = 2500, chains = 1, 
    control = list(adapt_delta = 0.99)
  )

pairs(fit, pars = c("sdo", "sdp", "lp__", "ar", "b"))

#-------------------------------------------------------------------------------
# plot stuff

p1 <- fit %>%
  spread_draws(R[year]) %>%
  median_qi() %>%
  ggplot(aes(x = year, y = R, ymin = .lower, ymax = .upper)) +
  geom_pointinterval()
p1
dat <- data.frame(year = 1:(n - k), Rtrue = R[(n_year - (n - (k + 1))):n_year])
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
dat <- data.frame(year = 1:n, Strue = S[(n_year - (n - 1)):n_year])
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
p3

p <- plot_grid(p3, p1, p2, ncol = 1)
p
# ggsave("plots/self-test-su-peterson.pdf", width = 8, height = 11)

# create a few posteriors for stock-recruit, i.e., ln(R/S) vs. S

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
  ylab(expression(Ln ~ frac(R, S))) +
  ggqfc::theme_qfc()

p4 <- p4 + geom_abline(
  intercept = devs$ar, slope = -devs$b, color = "black",
  linetype = 1, size = 0.5, alpha = .05
)

p4 <- p4 + geom_abline(
  intercept = ar_est, slope = -b_est, color = "white",
  linetype = 2, size = 1.5
)

p4 <- p4 + geom_abline(
  intercept = ar, slope = -b, color = "red",
  linetype = 1, size = 1.5
)

p4
