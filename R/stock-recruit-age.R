library(TMB)
library(tidyverse)
library(cowplot)
devtools::install_github("ChrisFishCahill/gg-qfc")
library(ggqfc)

# TODO: set up Req yet for Umsy + reference point calcs

# set leading parameters
Linf <- 1
vbk <- 0.2
alw <- 30
blw <- 3
lmat <- 0.6
ahv <- 3
sdv <- 0.7
Sad <- 0.8
Fmort <- 0.35
U <- 0.316
Ma <- -log(Sad)
ages <- 1:15
n_year <- 100
Ro <- 1
RecK <- 5 # Goodyear compensaiton ratio
sdwt <- 1.0 # standard deviation Rec
rhor <- 0 # autocorrelation for Rec
sdvt <- sqrt(0.1) # obs error in spawning stock size

# set up leading vectors
la <- Linf * (1 - exp(-vbk * ages))
wa <- alw * la^blw
fa <- wa - alw * lmat^blw
fa[which(fa < 0)] <- 0
vul <- 1 / (1 + exp(-1.7 * (ages - ahv) / sdv))
Snat <- Sad * (1 - exp(-1.5 * (ages - 0.8)))
# NOTE:
# just picked some reasonable values for the 1.7, 1.5 and 0.8
# but could be whatever

Surv0 <- SurvF <- rep(NA, length(ages))
Surv0[1] <- SurvF[1] <- 1
for (a in 2:length(ages)) {
  Surv0[a] <- Surv0[a - 1] * Snat[a - 1]
  SurvF[a] <- SurvF[a - 1] * Snat[a - 1] * (1 - U * vul[a - 1])
}

EPRo <- sum(fa * Surv0)
EPRf <- sum(fa * SurvF)
reca <- RecK / EPRo
recb <- -log(1 / RecK) / (Ro * EPRo)
sdc <- sqrt(1 - rhor^2) * sdwt
Bo <- Ro * EPRo
unfished <- fa * Surv0 / EPRo
fished <- fa * SurvF / EPRf
LorenzenS <- ((exp(vbk * (ages + 1)) - 1) / (exp(vbk * ages) - 1))^(-Ma / vbk)

# initialize recruits, eggs, nta matrix, etc
nta <- matrix(NA, nrow = n_year, ncol = length(ages))
eggs <- Ut <- rdev <- vulb <- ln_RS <- yield <- C <- rep(NA, n_year)
nta[1, ] <- Ro * Surv0
eggs[1] <- sum(nta[1, ] * fa)
vulb[1] <- sum(nta[1, ] * vul * wa)
yield[1] <- Ut[1] * sum(nta[1, ] * vul * wa)
vbo <- sum(nta[1, ] * wa * vul)

run_model <- function() {
  # draw devs & set up rec sequences
  rstd <- rnorm(n_year)
  rdev[1] <- sdc * rstd[1]
  for (t in 2:n_year) rdev[t] <- rhor * rdev[t - 1] + sdc * rstd[t]

  # set up exploitation history
  relU <- seq(from = 0, to = 1, by = 0.05)
  Ut[1:length(relU)] <- relU
  Ut[which(is.na(Ut))] <- 1
  Ut <- Ut * U
  C[1] <- Ut[1] * sum(nta[1, ] * vul)
  yield[1] <- Ut[1] * sum(nta[1, ] * vul * wa)

  # run the simulation
  for (t in 2:n_year) {
    for (a in 2:length(ages)) {
      nta[t, a] <- nta[t - 1, a - 1] * Snat[a - 1] * (1 - vul[a - 1] * Ut[t - 1])
      if (a == length(ages)) {
        nta[t, 1] <- reca * eggs[t - 1] * exp(-recb * eggs[t - 1] + rdev[t - 1])
        eggs[t] <- sum(nta[t, ] * fa)
        vulb[t] <- sum(nta[t, ] * vul * wa)
        yield[t] <- Ut[t] * sum(nta[t, ] * vul * wa)
        C[t] <- Ut[t] * sum(nta[t, ] * vul)
      }
    }
  }
  # true values
  r <- nta[2:n_year, 1]         # r from 2:n_year; recby in carl's code
  e <- eggs[1:(n_year-1)]       # eggs from 1:(n_year-1)
  ln_rs = log(r/e)              # true relationship

  # generate observed values
  vt <- rnorm(length(e), sdvt)  # iid sampling errors
  S <- e * exp(vt)              # Escapement w/ lognormal sampling error
  R <- S + C[1:(n_year - 1)]    # assume catch measured perfectly

  out <- tibble(
    # "true" stuff
    s = e, r = r, ln_rs = ln_rs,  
    yield = yield[1:(n_year - 1)],
    Ut = Ut[1:n_year - 1],
    vulb = vulb[1:(n_year - 1)], reca = rep(reca, n_year - 1),
    recb = rep(recb, n_year - 1), rdev = rdev[1:(n_year - 1)],
    Bo = rep(Bo, n_year - 1),
    sdc = rep(sdc, (n_year - 1)),  sdwt = rep(sdwt, (n_year - 1)), 
    sdvt = rep(sdvt, (n_year - 1)),
    rhor = rep(rhor, (n_year - 1)), 
    vt = vt, wt = rdev,
    year = 1:(n_year - 1), 
    C = C[1:(n_year - 1)],
    
    # corrupted stuff
    R = R, S = S, ln_RS = log(R/S)
  )
  out
}

# call the sim f(x) once to visualize
set.seed(3)
dat <- run_model()

p1 <- dat %>% ggplot(aes(x = year, y = yield)) +
  geom_point() +
  geom_line() +
  ggtitle("Yield") +
  theme_qfc() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- dat %>% ggplot(aes(x = year, y = vulb)) +
  geom_point() +
  geom_line() +
  ggtitle("Vulnerable biomass") +
  theme_qfc() +
  geom_hline(
    yintercept = vbo, linetype = 2,
    color = "steelblue", size = 0.7
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- dat %>% ggplot(aes(x = e, y = ln_rs)) +
  geom_point() +
  ylab("ln(r/s)") +
  ggtitle("ln(r/s) vs. s true") +
  theme_qfc() +
  theme(plot.title = element_text(hjust = 0.5))

p4 <- dat %>% ggplot(aes(x = S, y = ln_RS)) +
  geom_point() +
  ggtitle("ln(R/S) vs. S observed") +
  theme_qfc() +
  theme(plot.title = element_text(hjust = 0.5))

p5 <- dat %>% ggplot(aes(x = year, y = rdev)) +
  geom_point() +
  geom_line() +
  ggtitle("Recruitment anomalies") +
  theme_qfc() +
  theme(plot.title = element_text(hjust = 0.5))

p6 <- dat %>% ggplot(aes(x = year, y = Ut)) +
  geom_point() +
  geom_line() +
  ggtitle("Exploitation history") +
  theme_qfc() +
  theme(plot.title = element_text(hjust = 0.5))

p <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
p
ggsave("plots/sim-demonstration.pdf", width = 8, height = 10)

# demonstrate bias with dead simple linear regression
# note this is due to both errors in variables and time series bias
nsim <- 10000
ar_ests <- b_ests <- rep(NA, nsim)
set.seed(1)
for (i in 1:nsim) {
  dat <- run_model()
  fit <- lm(dat$ln_rs ~ dat$s) # ln(R/S) = ln(alpha) + beta*S
  ar_est <- fit$coefficients[1]
  b_est <- - fit$coefficients[2] # convert to -beta
  ar_ests[i] <- ar_est
  b_ests[i] <- b_est
}

# plot it
dat <- tibble(a_est = exp(ar_ests), b_est = b_ests)
a <- dat %>%
  ggplot(aes(x = a_est)) +
  geom_histogram(bins = 35) +
  geom_vline(
    xintercept = reca,
    color = "steelblue",
    linetype = 2, size = 1
  ) +
  xlab("Estimates of alpha vs. truth") +
  theme_qfc()

b <- dat %>%
  ggplot(aes(x = b_est)) +
  geom_histogram(bins = 35) +
  geom_vline(
    xintercept = recb,
    color = "steelblue",
    linetype = 2, size = 1
  ) +
  xlab("Estimates of beta vs. truth") +
  theme_qfc()

p <- plot_grid(a, b, ncol = 2)
p

ggsave("plots/reg-test.pdf", width = 7, height = 4)

#-----------------------
# TMB
# NOTE!!!
# Ths section a work in progress!

cppfile <- "src/sr_bias_simple.cpp"
compile(cppfile)
dyn.load(dynlib("src/sr_bias_simple"))

phi_start <- 0
parameters <- list(
  "ln_a" = log(1.6),
  "ln_b" = log(10.5),
  "ln_sdp" = log(0.05), "ln_sdm" = log(0.1),
  logit_phi = exp(phi_start) / (1 + exp(phi_start)),
  "ln_o" = rep(0, n_year),
  ln_ro = log(1)
)
random <- c("ln_o")

data <- list("R" = recby, "S" = S)
obj <- MakeADFun(data = data, parameters = parameters, random = random)

obj$fn(obj$par)
obj$gr(obj$par)

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)
report <- obj$report()
opt$SD <- sdreport(obj)
opt$SD

reca
a_ss[i] <- exp(unname(opt$par["log_a"]))
b_ss[i] <- exp(unname(opt$par["log_b"]))

plot(as.list(opt$SD, "Estimate")$wt)

par(mfrow = c(1, 2))
hist(a_ss,
  main = "state space ricker", xlab = "alpha est", breaks = 30,
  xlim = c(1, 5.5)
)
abline(v = reca, lwd = 2, col = "blue")

hist(b_ss, main = "state space ricker", xlab = "beta est", breaks = 30)
abline(v = recb, lwd = 2, col = "blue")
