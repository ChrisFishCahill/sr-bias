# simplest time series bias example
# cahill 29 Oct 2022
# set the true parameters
a <- 1   # actually ln(a) 
b <- 1
sdr <- 0.6
U <- 0.6
Req <- (a + log(1 - U)) / (b * (1 - U))

# set the required vectors
n_year <- 25
R <- S <- w <- ln_RS <- rep(NA, n_year)

# set up the TMB stuff: 
# compile and load the cpp
library(TMB)
cppfile <- "src/sr_bias_simple.cpp"
compile(cppfile)
dyn.load(dynlib("src/sr_bias_simple"))

# Build inputs
# data <- list("ln_RS" = ln_RS[-n_year], "S" = S[-n_year])
parameters <- list("log_a" = log(1), "log_sd_wt" = log(1), "log_sdm" = log(0.2), "log_b" = 1, "wt" = rep(0, n_year - 1))
random <- c("wt")
# use restricted maximum likelihood:
# random <- union(random, c("log_a", "log_b"))
# obj$fn(obj$par)
# obj$gr(obj$par)

# initialize first time step
R[1] <- Req
S[1] <- R[1] * (1 - U)

n_sims <- 10000
a_ests <- b_ests <- a_ss <- b_ss <- rep(NA, n_sims)
set.seed(1)
for (i in 1:n_sims) {
  norm_devs <- rnorm(n_year)
  wt <- sdr * norm_devs
  for (t in 2:n_year) {
    R[t] <- S[t - 1] * exp(a - b * S[t - 1] + wt[t - 1])
    S[t] <- R[t] * (1 - U)
    ln_RS[t - 1] <- log(R[t] / S[t - 1])
  }
  # normal ricker
  fit <- lm(ln_RS[-n_year] ~ S[-n_year])
  a_est <- fit$coefficients[1]
  b_est <- -fit$coefficients[2]
  a_ests[i] <- a_est
  b_ests[i] <- b_est
  
  # state-space ricker 
  data <- list("ln_RS" = ln_RS[-n_year], "S" = S[-n_year])
  obj <- MakeADFun(data = data, parameters = parameters, random = random)
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, control = list("trace" = 1))
  a_ss[i] = unname(opt$par["log_a"])
  b_ss[i] = unname(opt$par["log_b"])
  
  if(i%%100 == 0) cat(paste0("simulation \n i = ", i))
}

par(mfrow = c(2, 2))
hist(a_ests, main = "normal ricker", xlab = "ln alpha est", breaks = 30)
abline(v = a, lwd = 2, col = "blue")
hist(b_ests,
  main = "normal ricker", xlab = "beta est",
  xlim = c(0, 20), breaks = 70
)
abline(v = b, lwd = 2, col = "blue")

hist(a_ss, main = "state space ricker", xlab = "ln alpha est", breaks = 30)
abline(v = a, lwd = 2, col = "blue")

hist(b_ss, main = "state space ricker", xlab = "beta est", breaks = 30, 
     xlim = c(-5, 20))
abline(v = a, lwd = 2, col = "blue")
