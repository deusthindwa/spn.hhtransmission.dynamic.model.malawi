#written by Deus
#10/04/2026
#WAIFW household pneumococcal carriage transmission modelling in Malawi

#====================================================================

#stan model code

stan_code <- "

// ===========================================================================
// Pneumococcal WAIFW model — 5-group CTMC
// Group-specific clearance rates + FOI-dependent competition
//
// States:  1 = Susceptible,  2 = VT carrier,  3 = NVT carrier
// Groups:  1 = ychild, 2 = ochild, 3 = adult_hivneg,
//          4 = adult_hivS, 5 = adult_hivL
//
// KEY CHANGE vs previous version:
//   Competition rates are no longer fixed scalar parameters.
//   They are derived from the forces of infection at each transition:
//
//     sigma_VN(i,t) = lambda_NVT(i,t) * eps_V
//     sigma_NV(i,t) = lambda_VT(i,t)  * eps_N
//
//   eps_V = relative risk of NVT acquisition while currently carrying VT
//   eps_N = relative risk of VT  acquisition while currently carrying NVT
//
//   Biological consequence:
//     A VT carrier can only switch to NVT if there is active NVT exposure
//     (lambda_NVT > 0).  The switching rate scales with exposure pressure,
//     modulated by eps_V which captures competitive exclusion (eps < 1),
//     independence (eps = 1), or facilitation (eps > 1).
//     Because lambda_NVT and lambda_VT vary by individual and visit,
//     sigma_VN and sigma_NV are no longer constants — they are individual-
//     and visit-specific derived quantities inside the rate matrix Q.
// ===========================================================================

functions {

  /**
   * Build the 3x3 instantaneous rate matrix Q for individual n of group a.
   *
   * Competition rates are computed inside this function:
   *   sigma_VN = lNVT * epsV   (VT -> NVT switching)
   *   sigma_NV = lVT  * epsN   (NVT -> VT switching)
   *
   * This means sigma_VN and sigma_NV:
   *   (a) vary by individual and by transition (because lVT, lNVT do);
   *   (b) are zero whenever the competing type has zero exposure pressure;
   *   (c) scale with exposure to the competing type — capturing the
   *       mechanistic requirement that displacement requires encounter.
   *
   * @param lVT    total VT  FOI (S -> VT acquisition; used for sigma_NV)
   * @param lNVT   total NVT FOI (S -> NVT acquisition; used for sigma_VN)
   * @param muVT   VT  clearance rate of THIS individual's group (VT  -> S)
   * @param muNVT  NVT clearance rate of THIS individual's group (NVT -> S)
   * @param epsV   relative risk of NVT acquisition if currently VT carrier
   * @param epsN   relative risk of VT  acquisition if currently NVT carrier
   */
  matrix make_Q(real lVT,  real lNVT,
                real muVT, real muNVT,
                real epsV, real epsN) {

    // FOI-dependent competition rates for this individual at this transition
    real sigVN = lNVT * epsV;   // VT  -> NVT: NVT pressure x relative risk
    real sigNV = lVT  * epsN;   // NVT -> VT:  VT  pressure x relative risk

    matrix[3,3] Q;

    // Row 1: from Susceptible
    Q[1,1] = -(lVT + lNVT);
    Q[1,2] =   lVT;
    Q[1,3] =   lNVT;

    // Row 2: from VT carrier — clears at muVT[a]; switches to NVT at sigVN
    Q[2,1] =   muVT;
    Q[2,2] = -(muVT + sigVN);
    Q[2,3] =   sigVN;

    // Row 3: from NVT carrier — clears at muNVT[a]; switches to VT at sigNV
    Q[3,1] =   muNVT;
    Q[3,2] =   sigNV;
    Q[3,3] = -(muNVT + sigNV);

    return Q;
  }

}  // end functions


data {

  int<lower=1> N;       // number of individuals
  int<lower=1> N_AGE;   // number of age/HIV groups (= 5)

  array[N] int<lower=1,upper=3> s1;           // state at visit 1
  array[N] int<lower=1,upper=3> s2;           // state at visit 2
  array[N] int<lower=1,upper=3> s3;           // state at visit 3
  array[N] int<lower=1,upper=5> age_grp;      // group of each individual

  // HH carrier counts by age group [N x N_AGE x 2]
  // [individual, group, transition_index]
  array[N, N_AGE, 2] int<lower=0> n_VT;
  array[N, N_AGE, 2] int<lower=0> n_NVT;

  // Total other HH members [N x 2]
  array[N, 2] int<lower=0> n_HH_others;

  // Community prevalences [N_AGE x 2]  (age group, visit)
  matrix<lower=0,upper=1>[N_AGE, 2] comm_VT;
  matrix<lower=0,upper=1>[N_AGE, 2] comm_NVT;

}


parameters {

  // ---- WAIFW 5x5 household transmission matrices -------------------------
  // beta_H_VT[a,b]:  susceptible of group a acquires VT from infector of group b
  // beta_H_NVT[a,b]: same for NVT
  matrix<lower=0>[N_AGE, N_AGE] beta_H_VT;
  matrix<lower=0>[N_AGE, N_AGE] beta_H_NVT;

  // ---- Community force of infection [2 x N_AGE] --------------------------
  // lambda_C[1,a]: VT  community FOI for group a
  // lambda_C[2,a]: NVT community FOI for group a
  matrix<lower=0>[2, N_AGE] lambda_C;

  // ---- Group-specific clearance rates ------------------------------------
  // mu_VT[a]:  VT  clearance rate for carriers of group a
  // mu_NVT[a]: NVT clearance rate for carriers of group a
  vector<lower=0>[N_AGE] mu_VT;
  vector<lower=0>[N_AGE] mu_NVT;

  // ---- FOI-dependent competition parameters (relative risks) -------------
  //
  // eps_V: relative risk of acquiring NVT while currently carrying VT.
  //        Compared with the NVT acquisition rate of a fully susceptible
  //        individual (lambda_NVT), a VT carrier acquires NVT at rate
  //        lambda_NVT * eps_V.
  //        eps_V < 1 -> VT provides partial cross-protection against NVT
  //        eps_V = 1 -> no interaction (independent colonisation)
  //        eps_V > 1 -> VT carriage facilitates NVT acquisition
  //
  // eps_N: relative risk of acquiring VT while currently carrying NVT.
  //        Analogous interpretation.
  //
  // Prior: Exponential(2) => mean 0.5; encodes prior belief that competitive
  // exclusion (eps < 1) is more likely than independence or facilitation,
  // while permitting values above 1 if the data support them.
  real<lower=0> eps_V;
  real<lower=0> eps_N;

}


model {

  // ---- Priors ------------------------------------------------------------

  to_vector(beta_H_VT)  ~ normal(0.0, 1.0);   // half-normal (lower=0 constraint)
  to_vector(beta_H_NVT) ~ normal(0.0, 1.0);

  to_vector(lambda_C) ~ normal(0.0, 0.5);      // community FOI smaller than HH

  mu_VT  ~ normal(1.0, 0.5);    // ~1 clearance per inter-visit interval
  mu_NVT ~ normal(1.0, 0.5);

  // Relative-risk competition parameters
  // Exponential(2) => mean 0.5 (prior favour competitive exclusion, eps < 1)
  // allowing values > 1 (facilitation) if data support it
  eps_V ~ exponential(2.0);
  eps_N ~ exponential(2.0);

  // ---- Likelihood --------------------------------------------------------

  for (n in 1:N) {
    int a = age_grp[n];

    for (t in 1:2) {

      // ---- Forces of infection -------------------------------------------
      real denom = (n_HH_others[n,t] > 0) ? n_HH_others[n,t] : 1.0;

      // Initialise with community component
      real lVT  = lambda_C[1, a] * comm_VT[a,  t];
      real lNVT = lambda_C[2, a] * comm_NVT[a, t];

      // Accumulate household WAIFW contributions across all infector groups
      for (b in 1:N_AGE) {
        lVT  += beta_H_VT[a,  b] * n_VT[n,  b, t] / denom;
        lNVT += beta_H_NVT[a, b] * n_NVT[n, b, t] / denom;
      }

      // ---- Rate matrix with FOI-dependent competition --------------------
      // make_Q internally computes:
      //   sigma_VN(i,t) = lNVT * eps_V
      //   sigma_NV(i,t) = lVT  * eps_N
      matrix[3,3] Q = make_Q(
        lVT, lNVT,
        mu_VT[a], mu_NVT[a],
        eps_V, eps_N
      );

      // ---- Transition probabilities and log-likelihood -------------------
      matrix[3,3] P = matrix_exp(Q);

      int from_s = (t == 1) ? s1[n] : s2[n];
      int to_s   = (t == 1) ? s2[n] : s3[n];

      target += log(fmax(P[from_s, to_s], 1e-12));
    }
  }

}


generated quantities {

  // ---- Household reproduction numbers ------------------------------------
  // R_HH_VT[a,b]  = beta_H_VT[a,b]  / mu_VT[b]   (infector b's clearance)
  // R_HH_NVT[a,b] = beta_H_NVT[a,b] / mu_NVT[b]
  matrix[N_AGE, N_AGE] R_HH_VT;
  matrix[N_AGE, N_AGE] R_HH_NVT;

  // ---- Group-specific carriage durations ---------------------------------
  vector[N_AGE] dur_VT;    // 1 / mu_VT[a]
  vector[N_AGE] dur_NVT;   // 1 / mu_NVT[a]

  // ---- Competition summaries ---------------------------------------------
  //
  // relative_competition = eps_V / eps_N
  //   > 1: NVT displaces VT more readily (per unit of exposure pressure)
  //        than VT displaces NVT  =>  NVT competitive advantage
  //   < 1: VT displaces NVT more readily => VT competitive advantage
  //   = 1: symmetric competitive susceptibility
  //
  // Because sigma_VN(i,t) = lNVT(i,t) * eps_V and
  //         sigma_NV(i,t) = lVT(i,t)  * eps_N,
  // the realised ratio sigma_VN / sigma_NV also depends on the ratio
  // lNVT / lVT at each transition, which varies across individuals.
  // relative_competition isolates the competition parameters from exposure.
  real relative_competition;

  // Representative sigma values at the dataset-level mean FOI.
  // Computed as: mean(lVT) and mean(lNVT) averaged over all individual-
  // transition pairs, then multiplied by eps.
  // These scalars give an intuitive comparison with the old fixed-sigma model.
  real mean_sigma_VN;   // mean(lNVT) * eps_V
  real mean_sigma_NV;   // mean(lVT)  * eps_N

  // ---- Posterior predictive replication ----------------------------------
  array[N] int s2_rep;
  array[N] int s3_rep;

  // Compute R_HH, durations, and competition summaries
  relative_competition = eps_V / (eps_N + 1e-12);

  {
    real sum_lVT  = 0.0;
    real sum_lNVT = 0.0;
    int  cnt      = 0;

    for (n in 1:N) {
      int a = age_grp[n];

      for (t in 1:2) {
        real denom = (n_HH_others[n,t] > 0) ? n_HH_others[n,t] : 1.0;
        real lVT  = lambda_C[1, a] * comm_VT[a,  t];
        real lNVT = lambda_C[2, a] * comm_NVT[a, t];
        for (b in 1:N_AGE) {
          lVT  += beta_H_VT[a,  b] * n_VT[n,  b, t] / denom;
          lNVT += beta_H_NVT[a, b] * n_NVT[n, b, t] / denom;
        }
        sum_lVT  += lVT;
        sum_lNVT += lNVT;
        cnt += 1;
      }
    }

    real mean_lVT  = sum_lVT  / cnt;
    real mean_lNVT = sum_lNVT / cnt;

    mean_sigma_VN = mean_lNVT * eps_V;
    mean_sigma_NV = mean_lVT  * eps_N;
  }

  for (a in 1:N_AGE) {
    dur_VT[a]  = 1.0 / mu_VT[a];
    dur_NVT[a] = 1.0 / mu_NVT[a];
    for (b in 1:N_AGE) {
      R_HH_VT[a,b]  = beta_H_VT[a,b]  / mu_VT[b];
      R_HH_NVT[a,b] = beta_H_NVT[a,b] / mu_NVT[b];
    }
  }

  // Posterior predictive samples
  for (n in 1:N) {
    int a = age_grp[n];

    // t = 1 (visit 1 -> visit 2)
    {
      real denom = (n_HH_others[n,1] > 0) ? n_HH_others[n,1] : 1.0;
      real lVT  = lambda_C[1, a] * comm_VT[a, 1];
      real lNVT = lambda_C[2, a] * comm_NVT[a, 1];
      for (b in 1:N_AGE) {
        lVT  += beta_H_VT[a,  b] * n_VT[n,  b, 1] / denom;
        lNVT += beta_H_NVT[a, b] * n_NVT[n, b, 1] / denom;
      }
      matrix[3,3] P1 = matrix_exp(
        make_Q(lVT, lNVT, mu_VT[a], mu_NVT[a], eps_V, eps_N));
      s2_rep[n] = categorical_rng(to_vector(P1[s1[n]]));
    }

    // t = 2 (visit 2 -> visit 3)
    {
      real denom = (n_HH_others[n,2] > 0) ? n_HH_others[n,2] : 1.0;
      real lVT  = lambda_C[1, a] * comm_VT[a, 2];
      real lNVT = lambda_C[2, a] * comm_NVT[a, 2];
      for (b in 1:N_AGE) {
        lVT  += beta_H_VT[a,  b] * n_VT[n,  b, 2] / denom;
        lNVT += beta_H_NVT[a, b] * n_NVT[n, b, 2] / denom;
      }
      matrix[3,3] P2 = matrix_exp(
        make_Q(lVT, lNVT, mu_VT[a], mu_NVT[a], eps_V, eps_N));
      s3_rep[n] = categorical_rng(to_vector(P2[s2[n]]));
    }
  }

}
" 
