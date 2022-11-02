# simulate su and peterson
# cahill november 2022

k <- 2 # age at maturity
n_year <- 100
a <- 1.64 # ricker alpha
ar <- log(a)
b <- 1
sdp <- 0.01 # process error
sdo <- 0.3 # observation error
h <- runif(n_year - k, 0.01, 0.35) # harvest rate

E <- S <- rep(NA, n_year)
C <- R <- rep(NA, n_year - k) 

# Initialize S and C
S[1:k] <- ar / b # R' = S' = ln(a)/b = equilibrium SR

set.seed(1)
wt <- rnorm((n_year - 2), 0, sd = sdp)
for (t in 1:(n_year - k)) {
  R[t] <- a * S[t] * exp(-b * S[t] + wt[t]) # process noise 
  C[t] <- h[t] * R[t]
  S[t + k] <- R[t] - C[t]
}

# calculate observed values | process component of model
vt <- rnorm(n_year, 0, sdo)
E <- S * exp(vt) # obs. noise

plot(log(R / S[(k + 1):n_year]) ~ S[(k + 1):n_year],
     ylab = "ln(R/S) vs. S",
     xlab = "S"
)
plot(R ~ S[(k+1):n_year], "ylab" = "Recruits", xlab = "Stock")
summary(lm(log(R / S[(k + 1):n_year]) ~ S[(k + 1):n_year]))
plot(R ~ E[(k+1):n_year], "ylab" = "Recruits", xlab = "Escapement")

#-------------------------------------------------------------------------------
# TMB
#-------------------------------------------------------------------------------
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
