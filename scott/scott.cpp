#include <TMB.hpp>

template <class Type>
Type logitSelectivity(int a, Type mu, Type ups){
  Type tmp = Type(1) / (Type(1) + exp(Type(-1.0) * ups * (a - mu)));
  return tmp;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data inputs ---------------------------
  DATA_VECTOR(ages);
  DATA_VECTOR(paaC);      // prop at age female - commercial
  DATA_SCALAR(mu_M);      // M prior female- not log space
  DATA_SCALAR(sd_M);      // M prior female sd
  DATA_SCALAR(mu_F);      // F prior female - not log space
  DATA_SCALAR(sd_F);      // F prior female sd
  // DATA_SCALAR(mu_F50);    // sprF prior female - not log space
  // DATA_SCALAR(sd_F50);    // sprF prior female sd
  DATA_SCALAR(mu_mu);     // prior on selectivity S50
  DATA_SCALAR(sd_mu);     // prior on selectivity S50 sd
  DATA_SCALAR(mu_ups);    // prior on selectivity slope 
  DATA_SCALAR(sd_ups);    // prior on selectivity slope sd
  DATA_VECTOR(laa);       // length-at-age
  DATA_VECTOR(maa);       // maturrity-at-age
  DATA_VECTOR(waa);       // weight-at-age
  // DATA_SCALAR(target_spr);
  
  // Parameters -----------------------------
  PARAMETER(logM);		    // log nat mort rate female
  PARAMETER(logF); 		    // log female F commercial current
  PARAMETER(logmu);		    // log age full selectivity females, comm
  PARAMETER(logupsilon);	// log scale selectivity curve, comm
  PARAMETER(logsigR);
  // PARAMETER(logF50); 		   // log female F commercial current
  
  // Setup, procedures, init -------------
  
  int A = ages.size(); 	// ages a=1,...,A=30
  Type M = exp(logM);
  Type F = exp(logF);
  // Type F50 = exp(logF50);
  Type mu = exp(logmu);
  Type upsilon = exp(logupsilon);
  Type sigR = exp(logsigR);
  
  Type expnegM = exp(-M);
  vector<Type> Na(A); // Naa 
  vector<Type> Va(A); // unfished 
  vector<Type> Fa(A);
  // vector<Type> F50a(A);
  // vector<Type> SPRNa(A); // Naa 
  
  vector<Type> saC(A); // selectivity females - comm
  vector<Type> Ca(A); // catch females
  vector<Type> propC(A); // prop females - comm
  
  
  Type nll = 0.0; //init neg loglik
  // Type spr_nll = 0.0;
  Type tot_nll = 0.0; //init neg loglik
  
  // Priors --------------------------
  nll -=dnorm(M, mu_M, sd_M, true);
  nll -=dnorm(F, mu_F, sd_F, true);
  // nll -=dnorm(F50, mu_F50, sd_F50, true);
  nll -=dnorm(mu, mu_mu, sd_mu, true);
  nll -=dnorm(upsilon, mu_ups, sd_ups, true);
  
  // State dynamics ----------------------
  
  
  for(int a=0; a<A; a++){
    
    saC(a) = logitSelectivity(a+1, mu, upsilon);
    
  }
  
  Type maxsa = max(saC);
  for(int a=0; a<A; a++){
    saC(a) = saC(a) / maxsa;
  }
  
  for(int a=0; a<A; a++){
    Fa(a) = saC(a) * F;
  }
  
  // for(int a=0; a<A; a++){
    // F50a(a) = saC(a) * F50;
  // }
  
  Na(0) = 1000 ;
  Va(0) = 1000; 
  // SPRNa(0) = 1000;
  
  for(int a=1; a<(A-1); a++){
    Na(a) = Na(a-1) * exp(Type(-1.0) * (M + Fa(a-1))); 
    Va(a) = Va(a-1) * expnegM;
    // SPRNa(a) = SPRNa(a-1) * exp(Type(-1.0) * (M + F50a(a-1))); 
  }
  
  Na(A-1) = (Na(A-2) * exp(Type(-1.0) * (M + Fa(A-2)))) / (Type(1.0) - exp(Type(-1.0) * (M + Fa(A-1))));
  Va(A-1) = (Va(A-2) * expnegM) / (Type(1.0) - expnegM);
  // SPRNa(A-1) = (SPRNa(A-2) * exp(Type(-1.0) * (M + F50a(A-2)))) / (Type(1.0) - exp(Type(-1.0) * (M + F50a(A-1))));
  
  for(int a=0; a<A; a++){
    Ca(a) = Na(a) * (Type(1.0) - exp(Type(-1.0) * M - Fa(a))) * Fa(a) / (M + Fa(a));
  }
  
  Type maxf = sum(Ca);
  
  for(int a=0; a<A; a++){
    propC(a) = Ca(a) / maxf;
  }
  
  for(int a=0; a<A; a++){
    nll -=dnorm(propC(a), paaC(a), sigR, true);
    
  }
  
  
  Type unfished = 0;
  Type fished = 0;
  // Type SPRfished = 0;
  
  for(int a=0; a<A; a++){
    unfished += Va(a) * waa(a) * maa(a);
    fished += Na(a) * waa(a) * maa(a);
    // SPRfished += SPRNa(a) * waa(a) * maa(a);
  }
  
  Type spr = fished / unfished;
  
  // spr_nll = spr - target_spr;
  
  
  tot_nll += nll;
  // tot_nll += spr_nll;
  
  // Reports -------------------------------------
  
  REPORT(M);
  REPORT(F);
  REPORT(Fa);
  REPORT(saC);
  REPORT(propC);
  REPORT(mu);
  REPORT(upsilon);
  REPORT(Ca);
  REPORT(Na);
  REPORT(Va);
  // REPORT(F50);
  REPORT(spr);
  // REPORT(spr_nll);
  
  return tot_nll;
  
}
