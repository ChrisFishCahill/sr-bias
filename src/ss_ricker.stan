data {
  int<lower=0> k;
  int<lower=0> n_year;
  vector[n_year] E;
  vector[n_year-k] C; 
  vector[k] ar_prior; 
  vector[k] b_prior; 
  vector[k] sdp_prior; 
  vector[k] sdo_prior; 
  vector[k] So_prior; 
}
parameters {
  real ar; 
  real b; 
  real<lower=0> sdp;
  real<lower=0> sdo; 
  vector<lower=C>[n_year-k] R; 
  vector<lower=0>[k] So; 
}
transformed parameters {
 vector[n_year - k] mu; // log(R)
 vector[n_year] S;
 
 // initialize
 for(t in 1:n_year){S[t] = 0;}
 S[1] = So[1]; 
 S[2] = So[2]; 

 for(t in 1:(n_year-k)){
  mu[t] = ar + log(S[t]) - b*S[t];   
  S[t+k] = R[t] - C[t]; // R and C are t + k NOT t
 }
 
}
model {
  R ~ lognormal(mu, sdp);
  E ~ lognormal(log(S), sdo);
  // note, su and peterson use MVN for ar and b
  ar ~ normal(ar_prior[1], ar_prior[2]); 
  b ~ normal(b_prior[1], b_prior[2]); 
  sdp ~ normal(sdp_prior[1], sdp_prior[2]); 
  sdo ~ normal(sdo_prior[1], sdo_prior[2]); 
  So ~ lognormal(So_prior[1], So_prior[2]); 
}

