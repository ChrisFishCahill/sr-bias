#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(ln_RS);            // vector of R
  DATA_VECTOR(S);            // vector of S
  Type n_t = ln_RS.size(); 
  
  // Parameters
  PARAMETER(log_a);
  PARAMETER(log_b);
  PARAMETER(log_sdp);       // process sd
  PARAMETER(log_sdm);       // measurement sd 
  PARAMETER(rho);           // ar(1) parameter 
  PARAMETER_VECTOR(wt);     // latent state random coefficients
  PARAMETER(wt0);           // initial wt state
  
  // Objective function
  Type jnll = 0;
  
  // Probability of random coefficients--wt a vector of latent states
  jnll -= dnorm(wt(0), wt0, exp(log_sdp), true);  
  for( int t=1; t<n_t; t++){
    jnll -= dnorm( wt(t), rho*wt(t-1), exp(log_sdp), true); 
  }
  
  // Probability of data conditional on fixed and random effect values
  for( int t=0; t<n_t; t++){
    jnll -= dnorm( ln_RS(t), log_a - exp(log_b)*S(t) + wt(t), exp(log_sdm), true );
  }
  
  // Reporting
  Type sdp = exp(log_sdp);
  Type sdm = exp(log_sdm);
  REPORT(sdp);
  REPORT(sdm);
  
  return jnll;
}

