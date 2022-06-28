// urchsize
// This Stan code executes a size structure urchin population model 
//  fit to surey data of observed density and size structure of urchins 
//   
functions {
  // user function to create projection matrix
  matrix makemat(int Ns, vector G, vector S) {
     matrix[Ns,Ns] M; M = rep_matrix(0,Ns,Ns) ;
     M[1:(Ns-1),1:(Ns-1)] = add_diag(M[1:(Ns-1),1:(Ns-1)], S[1:(Ns-1)] .* (1 -  G[1:(Ns-1)])) ;
     M[2:Ns,1:(Ns-1)] = add_diag(M[2:Ns,1:(Ns-1)], S[1:(Ns-1)] .* G[1:(Ns-1)]) ;
     return M;
  }
}
// Section 1. Data inputs for model
data {
  int<lower=1> Nyrs ;            // Total years of counts (if 2011-2020 then 10)
  int<lower=1> YrB ;             // Year before temporal break (if 2013 then 3)
  int<lower=1> Nsz ;             // Number size classes (assumed to be 10)
  int<lower=1> Nszobs ;          // N years with size class data (if 2011-2020 then 10)
  int<lower=1> Sy[Nszobs] ;      // year index for size distribution observations
  vector[Nszobs] SzDst_N ;       // number of urchins measured for size distributions
  simplex[Nsz] Obs_Sz[Nszobs] ;  // Observed size distributions
  int<lower=1> Nsurvey ;         // N years with density surveys 
  int<lower=1> Dy[Nsurvey] ;     // year index for density surveys 
  vector<lower=0>[Nsurvey] Density ;// Observed density of urchins by year 
  vector[Nsz-1] Gr ;             // Growth transition probs for size classes
  simplex[Nsz] SZinit ;          // Initial size distribution
  real R0pri ;                   // prior for log mean Recruitment (R)
  real N0pri ;                   // prior for initial abundance
  vector[Nsz] Zeros ;            // vector of zeros (to make recruitment vector)
  real gamma_pri ;               // baseline log hazards (prior)
  real alpha ;                   // par1, logit fxn of Prob detection 
  vector[Nsz] sizes;             // vector of sizes for Prob detection fxn
}
// Section 2. The parameters to be estimated 
parameters {
  real<lower=0,upper=25> phi ;    // inverse scale param for gamma distribution
  real<lower=0,upper=50> theta[2] ;// params for precision of dirichlet dist. 
  real<lower=0,upper=25> sigR[2] ;// variation in log recruitment
  real<lower=0,upper=25> sigS[2] ;// variation in log hazards
  real<lower=0,upper=25> sigP[2]; // variation in logit detection prob
  vector[Nyrs] epsR ;             // random effect, recruitment by year
  vector[Nyrs] epsS ;             // random effect, survival (log hazards) by year
  vector[Nyrs] epsP ;             // random effect, detection prob (logit) by year
  real<lower=0> Beta ;            // par2, logit fxn of Prob detection  
  real<lower=-5,upper=5> gamma0 ; // baseline hazards (for survival)
  real R0 ;                       // baseline recruitment (log)
} 
// Section 3. Additional transformed parameters, including key model dynamics
transformed parameters {
  vector[Nyrs] R ;                // recruitment 
  vector[Nyrs] S;                 // survival rate
  vector[Nsz] P[Nyrs] ;           // detection rate
  vector[Nyrs] N ;                // true density, all years
  vector[Nsz] n[Nyrs] ;           // true population vector, all years
  vector[Nsz] d[Nyrs] ;           // expected detected urchins by year, size class
  vector[Nyrs] D ;                  // total detectable urchins, by year
  matrix[Nsz,Nsz] M0 ;            // initial matrix
  vector[Nsz] n0      ;           // initil pop vector
  vector[Nsz] lgtPD ;             // baseline logit detection fxn
  vector[Nsz] S0 ;                // baseline survival rates
  vector[Nsz] R0v ;               // nuisance recruitment vector
  real lampre_exp ;
  // Initial Prob detection vector
  for(i in 1:Nsz){
    lgtPD[i] = alpha - Beta * square(1 / sizes[i]) ;
  }  
  S0 =  rep_vector(exp(-exp(gamma0 - sigS[1] * mean(epsS[1:YrB]))),Nsz) ;   
  R0v = Zeros ;
  R0v[1] = exp(R0 + sigR[1] * mean(epsR[1:YrB]) ) ;
  n0 = N0pri * SZinit ;
  M0 = makemat(Nsz,Gr,S0) ; 
  for (t in 1:50){
    vector[Nsz] n0p = n0 ;
    n0 = M0 * (n0p + R0v) + .00001 ;
  }
  S[1] = exp(-exp(gamma0 - sigS[1] * epsS[1])) ;
  S0 = rep_vector(S[1],Nsz) ;
  R[1] = exp(R0 + sigR[1] * epsR[1] );
  R0v[1] = R[1] ;
  M0 = makemat(Nsz,Gr,S0) ; 
  n[1] = M0 * (n0 + R0v) + .00001 ;
  N[1] = sum(n[1]) ;
  P[1] = inv_logit(lgtPD + sigP[1] * epsP[1]) ;
  d[1] = n[1] .* P[1] ;
  D[1] = sum(d[1]) ;
  // Loop throuh years
  for (t in 2:Nyrs){
    matrix[Nsz,Nsz] Mt;             // matrix, year t
    vector[Nsz] Rt      ;           // recruitment vector, yr t
    vector[Nsz] St      ;           // survival vector, yr t
    real sg_S ;
    real sg_R ;
    real sg_P ;
    sg_S = t < YrB+1 ? sigS[1] : sigS[2];
    sg_R = t < YrB-1 ? sigR[1] : sigR[2];
    sg_P = t < YrB+1 ? sigP[1] : sigP[2];
    S[t] = exp(-exp(gamma0 - sg_S * epsS[t])) ;
    St = rep_vector(S[t],Nsz)  ; 
    R[t] = exp(R0 + sg_R * epsR[t]) ;
    Rt = Zeros ;
    Rt[1] = R[t] ;
    P[t] = inv_logit(lgtPD + sg_P * epsP[t]) ;
    Mt = makemat(Nsz,Gr,St) ;
    n[t] = Mt * (n[t-1] + Rt) + .00001 ;
    N[t] = sum(n[t]) ;
    d[t] = n[t] .* P[t] ;
    D[t] = sum(d[t]) ;
  }
}
// Section 4. Estimating model parameters (drawing from probability distributions)
model {
  // A) Observed nodes:
  // -Uchin size distrib
  for(i in 1:Nszobs){
    // Obs_Sz[t] ~ gamma(d[t] .* to_vector(phiM[t]) , phiM[t]) ; 
    Obs_Sz[i] ~ dirichlet( (d[Sy[i]] ./ D[Sy[i]]) * (theta[1]*pow(SzDst_N[i],theta[2]))) ;
  } 
  Density ~ gamma(D[Dy] * phi, phi ) ;
  // B) Prior distributions for model parameters:
  // Hierarchical random effects:  
  epsR ~ normal(0,1) ;
  epsS ~ normal(0,1) ;  
  epsP ~ normal(0,1) ;
  // Base parameter priors:
  sigR ~ cauchy(0,2.5) ;
  sigS ~ cauchy(0,2.5) ; 
  sigP ~ cauchy(0,2.5) ;
  phi ~ cauchy(0,2.5) ;
  theta ~ cauchy(0,2.5) ;
  Beta ~ cauchy(0,2.5) ;
  R0 ~ cauchy(0,2.5) ;
  gamma0 ~ cauchy(0,2.5) ;
}
// Section 5. Derived parameters and statistics 
 generated quantities {
  real Rmn_pre = mean(R[2:(YrB-2)]) ;
  real Rmn_post = mean(R[(YrB-1):(Nyrs-1)]) ;
  real logRdff = log(Rmn_post/Rmn_pre) ;
  real Smn_pre = mean(S[2:YrB]) ;
  real Smn_post = mean(S[(YrB+2):(Nyrs-1)]) ;
  real logSdff = log(Smn_post/Smn_pre) ;
  vector[Nyrs] Pmax ;
  real Pmn_pre ;  
  real Pmn_post ; 
  real logPdff ;
  vector[Nsz] SzDst[Nyrs] ;
  real ynew[Nsurvey]  ;          // New observations (out of sample) for ppc  
  vector[Nsurvey] P_resid;       // Pearson residuals, observed data
  vector[Nsurvey] P_resid_new;   // Pearson residuals, new data
  real Tstat ;                    // Test statistic, observed data (chi-2)
  real Tstat_new ;                // Test statistic, new data (chi-2)
  real ppp ;                      // posterior predictive P-value 
  for(t in 1:Nyrs){
    Pmax[t] = P[t][Nsz] ;
    SzDst[t] = d[t] ./ D[t] ;    
  }
  for(i in 1:Nsurvey){
    // residuals for observed and new data
    ynew[i] = gamma_rng(D[Dy[i]] * phi, phi) ;
    P_resid[i] = square(Density[i] - D[Dy[i]]) / (D[Dy[i]] / phi ) ;
    P_resid_new[i] = square(ynew[i] - D[Dy[i]]) / (D[Dy[i]] / phi) ; 
  }
  Pmn_pre = mean(Pmax[2:YrB]);
  Pmn_post  = mean(Pmax[(YrB+2):(Nyrs-1)]); 
  logPdff = log(Pmn_post/Pmn_pre) ;
  Tstat = sum(P_resid) ;
  Tstat_new = sum(P_resid_new) ;
  ppp = Tstat > Tstat_new ? 1 : 0;
}
