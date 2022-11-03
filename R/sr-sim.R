# simulate su and peterson
# cahill, punt november 2022

k <- 2 # age at maturity
n_year <- 100
a <- 1.6 # ricker alpha
ar <- log(a)
b <- 1
sdp <- 0.0 # process error sd 
sdo <- 0.1 # observation error sd

E <- S <- rep(NA, n_year) # Escapement, Stock
C <- R <- h <- rep(NA, n_year) # Catch, Recruits, finite harvest rate

# Initialize S, R, C 
set.seed(13)
h[1:k] = 0.1 # constant initial harvest rate
h[(k+1):n_year] <- runif(n_year-k, 0.01, 0.6) # harvest rate
wt <- rnorm(n_year-k, 0, sd = sdp) # process noise 

S[1:k] <- ar / b # S' = ln(a)/b = equilibrium S 
R[1:k] <- ar / b # R' = ln(a)/b = equilibrium R
C[1:k] <- h[1:k]*R[1:k] # equilibrium C given constant h = 0.1 

# sequentially generate new recruits, catch, and spawners 
for (t in 1:(n_year - k)) {
  R[t + k] <- a * S[t] * exp(-b * S[t] + wt[t]) # truth + process noise 
  C[t+k] <- h[t+k] * R[t+k]
  S[t + k] <- R[t+k] - C[t+k]
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
summary(lm(log(R / S) ~ S))


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
