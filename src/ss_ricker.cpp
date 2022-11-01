/*
           State-space Ricker stock-recruit model
              based on Su and Peterson (2012)
             1 November 2022 - Cahill QFC MSU
 
                          Maths: 
 y(t) = ln(E(t)); where E is Escapement
 x(t) = ln(S(t)); where S(t) is "true" number of spawners
 r(t) = ln(Rt); where R(t) is "true" number of recruits
 ar = ln(a); where a is Ricker alpha
 
 Observation equation: 
 y(t) = x(t) + v(t) for t = 1,...,n_year
 where T is total number of years
 and v(t) ~ Normal(0, sdvt) and thus represents
 lognormal observation/sampling error
 
 Process equation 1:
 r(t) = ar + x(t-k) - b*exp(x(t-k)) where k is age at maturity
 
 Process equation 2: 
 x(t) = ln(exp(r(t)) - C(t)) where C(t) is catch, assumed measured perfectly
 */

#include <TMB.hpp>

template <class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log = 0)
{
  Type logres;
  logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(k);               // age at maturity 
  DATA_VECTOR(E);                // escapements
  int n_year = E.size();         // total number of years
 
  PARAMETER(ar);                 // ln(a)
  PARAMETER(b);                  // ricker b
  PARAMETER_VECTOR(ln_Sinit);    // initialize the with Sinit vector of .size() = k
  PARAMETER(ln_sd_E);            // observation error
  PARAMETER(ln_sd_R);            // process error
  
  vector<Type> S(n_year);        // true spawners
  vector<Type> R(n_year);        // true recruits
  vector<Type> mu(n_year);       // predictions
  
  S.setZero(); R.setZero(); mu.setZero();  
  
  // initialize state variables
  for(int t = 0; t < k; t ++){
    S(t) = exp(ln_Sinit(t));
  }
  
  Type jnll = 0;
  for(int t = 0; t < n_year - k; t++){
    mu(t) = ar + log(S(t)) - b*S(t); 
    S(t + k) = R(t) - E(t);                        // note R and E represent R(t+k) and E(t+k)
    jnll -= dlnorm(R(t), mu(t), exp(ln_sd_R));     // process error in recruitment
  }

  for(int t = 0; t < k; t ++){
    jnll -= dlnorm(E(t), log(S(t)), exp(ln_sd_E)); // observation error
  }
  
  return jnll;
}
