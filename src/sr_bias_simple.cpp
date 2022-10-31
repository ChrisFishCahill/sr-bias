#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(R); 
  DATA_VECTOR(S);            
  Type n_t = R.size(); 
  
  PARAMETER(ln_a);
  PARAMETER(ln_b);
  PARAMETER(ln_sdp);                // process sd
  PARAMETER(ln_sdm);                // measurement sd 
  PARAMETER(logit_phi);             // ar(1) parameter 
  PARAMETER_VECTOR(ln_o);           // latent state random coefficients
  PARAMETER(ln_ro); 
  Type phi = invlogit(logit_phi); 
  
  Type jnll = 0;
  // probability of random coefficients
  jnll -= dnorm(ln_o(0), ln_ro, exp(ln_sdp), true);  
  for( int t=1; t<n_t; t++){
    jnll -= dnorm( ln_o(t), log(S(t)) + ln_a - exp(ln_b)*S(t) + phi*ln_o(t-1), exp(ln_sdp), true); 
  }

  // probability of data conditional on fixed and random effect values
  for( int t=0; t<n_t; t++){
    jnll -= dnorm( log(R(t)), ln_o(t), exp(ln_sdm), true );
  }
  
  // reporting
  Type sdp = exp(ln_sdp);
  Type sdm = exp(ln_sdm);
  REPORT(sdp);
  REPORT(sdm);
  REPORT(phi); 
  
  return jnll;
}

