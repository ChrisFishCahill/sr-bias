#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(ln_RS);        // vector of ln(R/S)
  DATA_VECTOR(S);            // vector of S
  Type n_t = ln_RS.size(); 
  
  // Parameters
  PARAMETER(log_a);
  PARAMETER(log_sdwt);      // process sd
  PARAMETER(log_sdm);       // measurement sd 
  PARAMETER(log_b);
  PARAMETER_VECTOR(wt);     // latent state random coefficients
  
  // Objective function
  Type jnll = 0;
  
  // Probability of random coefficients--wt a vector of latent states
  for( int t=0; t<n_t; t++){
    jnll -= dnorm( wt(t), Type(0.0), exp(log_sdwt), true);
  }
  
  // Probability of data conditional on fixed and random effect values
  for( int t=0; t<n_t; t++){
    jnll -= dnorm( ln_RS(t), log_a - exp(log_b)*S(t) + wt(t), exp(log_sdm), true );
  }
  
  // Reporting
  Type sdwt = exp(log_sdwt);
  Type sdm = exp(log_sdm);
  REPORT(wt);
  REPORT(sdwt);
  REPORT(sdm);
  
  return jnll;
}

