model
{
  # Observation equation
  for (t in 1:T)
  {
    E[t] ∼ dlnorm(mu.logE[t], prec.E); mu.logE[t] <- log(S[t]);
  }
  # Prior distributions for the state variables
  for (t in 1:k)
  {
    S[t] <- S1[t]; S1[t] ∼ dlnorm(0.0, 1.0E−3)
  }
  for (t in 1:(T −k))
  {
    mu[t] <- b[1] + log(S[t]) − b[2] * S[t];
    R[t] ∼ dlnorm(mu[t], prec.R)I(C[t], 1.0E+9); # R and C represent Rt+k and Ct+k, respectively
    S[t + k] <- R[t] - C[t];
  }
  # Prior distributions for model parameters
  # a bivariate normal prior for a r and beta
  b[1:2] ∼ dmnorm(mu.b[], prec[,])
  ar <- b[1]; alpha <- exp(b[1]);
  beta <- b[2]
  # Priors for standard deviations
  tau  ∼ dunif(0, 1.0E+9); tau2 <- tau * tau; prec.R <- 1/tau2
  sig ∼ dunif(0, 1.0E+9); sig2 <- sig * sig; prec.E <-  1/sig2
  # Management-oriented quantities
  SMSY <- hMSY/beta
  hMSY <- 0.5*ar - 0.07*ar *ar
}