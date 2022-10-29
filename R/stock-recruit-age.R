# leading parameters
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

# leading vectors
la <- Linf * (1 - exp(-vbk * ages))
wa <- alw * la^blw
fa <- wa - alw * lmat^blw
fa[which(fa < 0)] <- 0
vul <- 1 / (1 + exp(-1.7 * (ages - ahv) / sdv))
Snat <- Sad * (1 - exp(-1.5 * (ages - 0.8)))
# where do the 1.7, 1.5 and 0.8 come from?
Surv0 <- SurvF <- rep(NA, length(ages))
Surv0[1] <- SurvF[1] <- 1
for (a in 2:length(ages)) {
  Surv0[a] <- Surv0[a - 1] * Snat[a - 1]
  SurvF[a] <- SurvF[a - 1] * Snat[a - 1] * (1 - U * vul[a - 1])
}

Ro <- 1
RecK <- 5
EPRo <- sum(fa * Surv0)
EPRf <- sum(fa * SurvF)
reca <- RecK / EPRo
recb <- -log(1 / RecK) / (Ro * EPRo)
sdr <- 0.6
rhor <- 0
sdc <- sqrt(1 - rhor^2) * sdr
Bo <- Ro * EPRo
# Req <- ? =IF(B35>0,MAX(0,-LN(1/(B34*EPRf*EXP(0.5*sdr^2)))/(B35*EPRf)),0)
unfished <- fa * Surv0 / EPRo
fished <- fa * SurvF / EPRf
LorenzenS <- ((exp(vbk * (ages + 1)) - 1) / (exp(vbk * ages) - 1))^(-Ma / vbk)

# initialize recruits, eggs, na matrix, etc
nta <- matrix(NA, nrow = n_year, ncol = length(ages))
r <- eggs <- Ut <- rdev <- vulb <- ln_RS <- yield <- rep(NA, n_year)
nta[1, ] <- Ro * Surv0
eggs[1] <- sum(nta[1, ] * fa)
vulb[1] <- sum(nta[1, ] * vul * wa)
yield[1] <- Ut[1] * sum(nta[1, ] * vul * wa)

# pick some deviates
set.seed(1)
rstd <- rnorm(n_year)
rdev[1] <- sdc * rstd[1]
for (t in 2:n_year) rdev[t] <- rhor * rdev[t - 1] + sdc * rstd[t]
relU <- seq(from = 0, to = 1, by = 0.05)
Ut[1:length(relU)] <- relU
Ut[which(is.na(Ut))] <- 1
Ut <- Ut * U

# run the model
for (t in 2:n_year) {
  for (a in 2:length(ages)) {
    nta[t, a] <- nta[t - 1, a - 1] * Snat[a - 1] * (1 - vul[a - 1] * Ut[t - 1])
    if (a == length(ages)) {
      nta[t, 1] <- reca * eggs[t - 1] * exp(-recb * eggs[t - 1] + rdev[t - 1])
      eggs[t] <- sum(nta[t, ] * fa)
      vulb[t] <- sum(nta[t, ] * vul * wa)
      yield[t] <- Ut[t] * sum(nta[t, ] * vul * wa)
    }
  }
}

S <- nta[1:(n_year - 1), 1]
ln_RS <- eggs[1:(n_year - 1)] / S

par(mfrow = c(1, 2))
plot(yield, type = "b")
plot(vulb, type = "b")
plot(ln_RS ~ S)

#-----------------------
library(TMB)

cppfile <- "src/sr_bias_simple.cpp"
compile(cppfile)
dyn.load(dynlib("src/sr_bias_simple"))
parameters <- list(
  "log_a" = log(1), "log_sdwt" = log(1),
  "log_sdm" = log(0.2),
  "log_b" = 1, "wt" = rep(0, n_year-1)
)
random <- c("wt")

data <- list("ln_RS" = ln_RS, "S" = S)
obj <- MakeADFun(data = data, parameters = parameters, random = random)
obj$fn(obj$par)
obj$gr(obj$par)

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, control = list("trace" = 1))
exp(unname(opt$par["log_a"]))
