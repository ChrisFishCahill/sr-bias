# simulate su and peterson
# cahill november 2022


# Su and peterman escapement data example
# Ep = c(1.062,
#       1.179, 1.385, 1.213, 1.176, 1.385, 1.141, 0.944, 1.087, 1.119,
#       0.966, 0.918, 0.793, 0.775, 0.645, 0.661, 0.548, 0.703, 0.807,
#       0.792, 0.9, 1.007, 1.222, 1.269, 1.001, 1.108, 1.032, 1.099)
#
#
# C = c(0.165, 0.23, 0.272, 0.34, 0.385, 0.409, 0.508, 0.596,
#        0.609, 0.669, 0.674, 0.779, 0.708, 0.864, 0.611, 0.658, 0.631,
#        0.583, 0.548, 0.534, 0.556, 0.513, 0.332, 0.319, 0.235, 0.2)

k <- 2 # age at maturity
n_year <- 100
a <- 1.64 # ricker alpha
ar <- log(a)
b <- 1

# make some fake catches -- taken from su and peterson
Sinit <- rep(1, k)
E <- S <- rep(NA, n_year)
C <- mu <- R <- rep(NA, n_year - k) # mu is log space predictor of R
sdp <- 0.03 # process error
h <- runif(n_year - k, 0.01, 0.35) # harvest rate

# Initialize S and C
S[1:k] <- ar / b # R' = S' = ln(a)/b = equilibrium SR

set.seed(1)
wt <- rnorm((n_year - 2), 0, sd = sdp)

for (t in 1:(n_year - k)) {
  R[t] <- a * S[t] * exp(-b * S[t] + wt[t])
  C[t] <- h[t] * R[t]
  S[t + k] <- R[t] - C[t]
}

plot(log(R / S[(k + 1):n_year]) ~ S[(k + 1):n_year],
  ylab = "ln(R/S) vs. S",
  xlab = "S"
)
# plot(R ~ S[(k+1):n_year], "ylab" = "Recruits", xlab = "Stock")
summary(lm(log(R / S[(k + 1):n_year]) ~ S[(k + 1):n_year]))

# calculate observed values based on values from process component of model
sdo <- 0.1 # observation error
vt <- rnorm(n_year, 0, sdo)
E <- S * exp(vt) # obs. noise

# TMB
library(TMB)

data <- list(
  "k" = k,
  "E" = E,
  "C" = C
)
par <- list(
  "ar" = ar,
  "b" = b,
  "ln_sdo" = log(sdo),
  "ln_sdp" = log(sdp),
  "R" = R
)

cppfile <- "src/ss_ricker.cpp"
compile(cppfile)
dyn.load(dynlib("src/ss_ricker"))

obj <- MakeADFun(data = data, parameters = par)

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

plot(obj$report()$`R` ~ R)
