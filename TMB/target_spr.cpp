#include <TMB.hpp>

// SPR model with single age composition

template <class Type> Type square(Type x){return x*x;} // user defined square() function

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data inputs ---------------------------
  DATA_VECTOR(ages);
  DATA_VECTOR(paaC);       // prop at age female - commercial
  DATA_SCALAR(mu_M);       // M prior female- not log space
  DATA_SCALAR(sd_M);       // M prior female sd
  DATA_SCALAR(mu_F);       // F prior female - not log space
  DATA_SCALAR(sd_F);       // F prior female sd
  DATA_SCALAR(mu_Fspr);    // sprF prior female - not log space
  DATA_SCALAR(sd_Fspr);    // sprF prior female sd
  DATA_SCALAR(mu_mu);      // prior on selectivity S50
  DATA_SCALAR(sd_mu);      // prior on selectivity S50 sd
  DATA_SCALAR(mu_ups);     // prior on selectivity slope 
  DATA_SCALAR(sd_ups);     // prior on selectivity slope sd
  DATA_VECTOR(laa);        // length-at-age
  DATA_VECTOR(maa);        // maturity-at-age
  DATA_VECTOR(waa);        // weight-at-age
  DATA_SCALAR(target_spr); // target SPR (e.g. 0.5)
  
  // Parameters -----------------------------
  PARAMETER(logM);		     // log nat mort rate female
  PARAMETER(logF); 		     // log female F commercial current
  PARAMETER(logmu);		     // log age at 50% selectivity females
  PARAMETER(logupsilon);	 // log selectivity slope
  PARAMETER(logsigR);      // log variability in comps
  PARAMETER(logFspr); 	   // log fishing mortality rate at target SPR
  
  // Setup, procedures, init -------------
  
  int A = ages.size(); 	   // number of ages
  Type M = exp(logM);
  Type F = exp(logF);
  Type Fspr = exp(logFspr);
  Type mu = exp(logmu);
  Type upsilon = exp(logupsilon);
  Type sigR = exp(logsigR);

  Type expnegM = exp(-M); // survival
  vector<Type> Na(A); // Naa 
  vector<Type> Va(A); // unfished 
  vector<Type> Fa(A);
  vector<Type> Fspra(A);
  vector<Type> SPRNa(A); // Naa
  
  vector<Type> saC(A); // selectivity females - comm
  vector<Type> Ca(A); // catch females
  vector<Type> propC(A); // prop females - comm
  
  
  // Likelihood components
  Type priors = 0.0; 
  Type comp_nll = 0.0;
  Type spr_nll = 0.0;
  Type tot_nll = 0.0; 
  
  // Priors -----------------------------
  
  priors -=dnorm(M, mu_M, sd_M, true);
  priors -=dnorm(F, mu_F, sd_F, true);
  priors -=dnorm(Fspr, mu_Fspr, sd_Fspr, true);
  priors -=dnorm(mu, mu_mu, sd_mu, true);
  priors -=dnorm(upsilon, mu_ups, sd_ups, true);
  
  // State dynamics ----------------------
  
  // Logistic selectivity-at-age
  for(int a=0; a<A; a++){
    saC(a) = Type(1) / (Type(1) + exp(Type(-1.0) * upsilon * (a - mu)));
  }

  // Fully-selected fishing mortality at age
  for(int a=0; a<A; a++){
    Fa(a) = saC(a) * F;
  }
  
  // TODO what's this?
  for(int a=0; a<A; a++){
  Fspra(a) = saC(a) * Fspr;
  }
  
  // Initialize values
  Na(0) = 1000;
  Va(0) = 1000; 
  SPRNa(0) = 1000;
  
  for(int a=1; a<(A-1); a++){
    // Fished numbers-at-age
    Na(a) = Na(a-1) * exp(Type(-1.0) * (M + Fa(a-1))); 
    // Virgin/unfished numbers-at-age
    Va(a) = Va(a-1) * expnegM;
    // Numbers-at-age at target SPR
    SPRNa(a) = SPRNa(a-1) * exp(Type(-1.0) * (M + Fspra(a-1)));
  }
  
  // Plus group dynamics
  Na(A-1) = (Na(A-2) * exp(Type(-1.0) * (M + Fa(A-2)))) / (Type(1.0) - exp(Type(-1.0) * (M + Fa(A-1))));
  Va(A-1) = (Va(A-2) * expnegM) / (Type(1.0) - expnegM);
  SPRNa(A-1) = (SPRNa(A-2) * exp(Type(-1.0) * (M + Fspra(A-2)))) / (Type(1.0) - exp(Type(-1.0) * (M + Fspra(A-1))));
  
  // Baranov catch equation to get catch in numbers-at-age
  for(int a=0; a<A; a++){
    Ca(a) = Na(a) * (Type(1.0) - exp(Type(-1.0) * M - Fa(a))) * Fa(a) / (M + Fa(a));
  }
  
  // Predicted catch-at-age composition (turn into a proportion)
  Type sumC = sum(Ca);
  for(int a=0; a<A; a++){
    propC(a) = Ca(a) / sumC;
  }
  
  // Likelihood component for age comp data TODO: not "normal" to use a normal
  // distribution for comp data - use multinomial instead?
  for(int a=0; a<A; a++){
    comp_nll -= dnorm(propC(a), paaC(a), sigR, true);
  }
  
  Type unfished = 0.0; // Initialize
  Type fished = 0.0;
  Type SPRfished = 0.0;
  
  // Spawning biomass potential (mutliply numbers-at-age by weight-at-age and
  // maturity-at-age)
  for(int a=0; a<A; a++){
    unfished += Va(a) * waa(a) * maa(a);
    fished += Na(a) * waa(a) * maa(a);
    SPRfished += SPRNa(a) * waa(a) * maa(a);
  }
  
  Type catch_target_spr = 0.0;
  // Baranov catch equation to get catch in numbers-at-age
  for(int a=0; a<A; a++){
    catch_target_spr += waa(a) * SPRNa(a) * (Type(1.0) - exp(Type(-1.0) * M - Fspra(a))) * Fspra(a) / (M + Fspra(a));
  }
  
  // spawning potential ratio
  Type spr = fished / unfished;
  
  // SPR penalty (if spr_target = 0.5, then minimizes ratio of
  // SPRfished/unfished to be 0.5)
  spr_nll = 100 * square(SPRfished / unfished - target_spr);
  
  // Add likelihood components
  tot_nll = priors + comp_nll + spr_nll;

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
  REPORT(Fspr);
  REPORT(spr);
  REPORT(catch_target_spr);
  REPORT(priors);
  REPORT(comp_nll);
  REPORT(spr_nll);
  REPORT(tot_nll);
  
  return tot_nll;
  
}
