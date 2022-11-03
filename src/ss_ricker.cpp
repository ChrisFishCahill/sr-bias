/*
   State-space Ricker stock-recruit model
      based on Su and Peterson (2012)
            Cahill QFC MSU
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
  DATA_INTEGER(k);                                        // age at maturity 
  DATA_VECTOR(E);                                         // escapements - observed with error, 1:n_year
  DATA_VECTOR(C);                                         // catches - assumed perfectly known
  int n_year = E.size();                                  // total number of years

  PARAMETER(ar);                                          // ln(a)
  PARAMETER(b);                                           // ln(b)
  PARAMETER(ln_sdo);                                      // observation error
  PARAMETER(ln_sdp);                                      // process error
  PARAMETER_VECTOR(R);                                    // recruits from (k+1):n_year
  PARAMETER_VECTOR(ln_So);                                // initial S(0:(k-1))
  // containers                         
  vector<Type> S(n_year);                                 // spawners
  vector<Type> mu(n_year-k);                              // ln(predicted recruits) from (k+1):n_year
  S.setZero(); mu.setZero();     
  
  Type sdo = exp(ln_sdo);    
  Type sdp = exp(ln_sdp);    
  Type arec = exp(ar);    

  // initialize 
  for(int t = 0; t < k; t ++){ S(t) = exp(ln_So(t));}    
  
  Type jnll = 0;
  
  for(int t = 0; t < n_year - k; t++){ 
    mu(t) = ar + log(S(t)) - b*S(t);     
    S(t + k) = R(t) - C(t);                               // R and C are t+k NOT t
    jnll -= dlnorm(R(t), mu(t), exp(ln_sdp), true);       // process error
  }   
   
  for(int t = 0; t < n_year; t ++){   
    jnll -= dlnorm(E(t), log(S(t)), exp(ln_sdo), true);   // observation error
  }
  
  REPORT(S); 
  REPORT(C); 
  REPORT(mu); 
  REPORT(E); 
  REPORT(sdp); 
  REPORT(sdo); 
  REPORT(arec); 

  return jnll;
}
