# simulate su and peterson
# cahill, punt november 2022

k <- 2 # age at maturity
n_year <- 100
ar <- 1 # ln( ricker alpha)
a <- exp(ar) 
b <- 1
sdp <- 0.05 # process error sd
sdo <- 0.1 # observation error sd

E <- S <- rep(NA, n_year) # Escapement, Stock
C <- R <- h <- rep(NA, n_year) # Catch, Recruits, finite harvest rate

# Initialize S, R, C
Ut <- rep(NA, n_year)
U = 0.45
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
#plot(S)
#plot(R)
#summary(lm(log(R / S) ~ S))
confint(lm(log(R / S) ~ S))

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
  "br" = log(b),
  "ln_sdo" = log(sdo),
  "ln_sdp" = log(sdp),
  "ln_R" = log(jitter(R[(k+1):n_year])),
  "ln_So" = c(log(S[1]), log(S[2]))
)

cppfile <- "src/ss_ricker.cpp"
compile(cppfile)
dyn.load(dynlib("src/ss_ricker"))
random = "ln_R"
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


