#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(ln_RS);            
  DATA_VECTOR(S);            
  Type n_t = ln_RS.size(); 
  
  PARAMETER(ln_a);
  PARAMETER(ln_b);
  PARAMETER(ln_sdp);                // process sd
  PARAMETER(ln_sdm);                // measurement sd 
  PARAMETER(logit_rho);             // ar(1) parameter 
  PARAMETER_VECTOR(wt);             // latent state random coefficients
  // PARAMETER(wto); 
  Type rho = invlogit(logit_rho); 
  
  Type jnll = 0;
  
  // probability of random coefficients--wt is a vector of latent states
  jnll -= dnorm(wt(0), ln_a - exp(ln_b)*S(0), exp(ln_sdp), true);  
  for( int t=1; t<n_t-1; t++){
    jnll -= dnorm( wt(t), rho*wt(t-1), exp(ln_sdp)/(1-pow(rho,2)), true); 
  }

  // probability of data conditional on fixed and random effect values
  for( int t=1; t<n_t; t++){
    jnll -= dnorm( ln_RS(t), ln_a - exp(ln_b)*S(t) + wt(t), exp(ln_sdm), true );
  }
  
  // reporting
  Type sdp = exp(ln_sdp);
  Type sdm = exp(ln_sdm);
  REPORT(sdp);
  REPORT(sdm);
  REPORT(rho); 
  
  return jnll;
}

