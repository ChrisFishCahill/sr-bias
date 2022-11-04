data {
  int<lower=0> k;
  int<lower=0> n_year;
  vector[n_year] E;
  vector[n_year-k] C; 
  vector[k] ar_prior; 
  vector[k] b_prior; 
  real sdp_prior; 
  real sdo_prior; 
}
parameters {
  real ar; 
  real b; 
  real<lower=0> sdp;
  real<lower=0> sdo; 
  vector<lower=C>[n_year-k] R; 
}
transformed parameters {
 vector[n_year - k] mu; 
 vector[n_year] S;
 
 // initialize
 for(t in 1:n_year){S[t] = ar/b;}

 for(t in 1:(n_year-k)){
  mu[t] = ar + log(S[t]) - b*S[t];   
  S[t+k] = R[t] - C[t]; // R and C are t + k NOT t
 }
}
model {
  // priors 
  // note, su and peterson use MVN for ar and b
  ar ~ normal(ar_prior[1], ar_prior[2]); 
  b ~ normal(b_prior[1], b_prior[2]); 
  sdp ~ exponential(sdp_prior); 
  sdo ~ exponential(sdo_prior); 

  // likelihoods
  R ~ lognormal(mu, sdp);
  E ~ lognormal(log(S), sdo);
}

