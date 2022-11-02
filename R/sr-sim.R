# simulate su and peterson
# xt = ln(St)
# rt = ln(Rt)
# ar = ln(alpha)
# y = xt + vt
# rt = ar + x(t-k) - b*x(t-k) + w(t-k)
# xt = ln(exp(rt) - C)
# need to initialize Sinit 1,2

k <- 2 # age at maturity
n_year = 28
ar = 1.1
a = exp(ar)
b = 3

# make some fake catches -- taken from su and peterson
C <- c(
  0.165, 0.23, 0.272, 0.34, 0.385, 0.409, 0.508, 0.596,
  0.609, 0.669, 0.674, 0.779, 0.708, 0.864, 0.611, 0.658, 0.631,
  0.583, 0.548, 0.534, 0.556, 0.513, 0.332, 0.319, 0.235, 0.2
)

Sinit = rep(2, k)
S = rep(NA, n_year)
mu <- R <- E <- rep(NA, n_year-k) # mu is log space predictor of R
sdp = 0.005
sdo = 0.01

# set up mu, R, and S
S[1:k] = Sinit # initialize
for(t in 1:(n_year - k)){
  mu[t] = ar + log(S[t]) - b*S[t]
  R[t] = rlnorm(1, exp(mu[t]), sdp) # process noise
  S[t + k] = R[t] - C[t]
}

# calculate escapements based on true value 
# S[t] and lognormal observation noise:
for(t in 1:n_year){
  E[t] = rlnorm(1, exp(log(S[t])), sdo) # obs. noise
}

# TMB
library(TMB)
data <- list(
  "k" = 2, 
  "E" = c(
    1.062,
    1.179, 1.385, 1.213, 1.176, 1.385, 1.141, 0.944, 1.087, 1.119,
    0.966, 0.918, 0.793, 0.775, 0.645, 0.661, 0.548, 0.703, 0.807,
    0.792, 0.9, 1.007, 1.222, 1.269, 1.001, 1.108, 1.032, 1.099
  ),
  "C" = c(
    0.165, 0.23, 0.272, 0.34, 0.385, 0.409, 0.508, 0.596,
    0.609, 0.669, 0.674, 0.779, 0.708, 0.864, 0.611, 0.658, 0.631,
    0.583, 0.548, 0.534, 0.556, 0.513, 0.332, 0.319, 0.235, 0.2
  )
)
par <- list(
  "ar" = log(1.4),
  "b" = 0.5,
  "ln_Sinit" = rep(log(1),2),
  "ln_sdp" = log(0.005),
  "ln_sdo" = log(0.01), 
  "ln_R" = log(R)
  )

cppfile <- "src/ss_ricker.cpp"
compile(cppfile)
dyn.load(dynlib("src/ss_ricker"))

obj <- MakeADFun(data = data, parameters = par)

obj$fn(obj$par)
obj$gr(obj$par)
obj$report()

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)



