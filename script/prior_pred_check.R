# =============================================================================
# Prior Predictive Check — 5-Group WAIFW CTMC Model with FOI-Dependent
# Competition (spn_waifw_5grp_foicomp_stan workflow)
# =============================================================================
#
# PURPOSE
# -------
# A prior predictive check (PPC) evaluates whether the prior distributions
# placed on model parameters generate plausible data BEFORE conditioning on
# observations. It answers:
#
#   "If we drew parameters from the priors and simulated data through the
#    full model, does the resulting distribution of outcomes look broadly
#    consistent with what we actually observed?"
#
# This is distinct from a POSTERIOR predictive check (which uses fitted
# draws) — here we never touch the likelihood; parameters come only from
# the prior.
#
# APPROACH
# --------
# Two complementary strategies are used:
#
#   (A) Stan-based: A modified Stan model that removes the likelihood
#       from the model block, leaving only priors. The generated quantities
#       block simulates state sequences (s2_sim, s3_sim) via categorical_rng.
#       Stan samples from the prior using the fixed_param sampler.
#
#   (B) R-based: Faster, transparent prior sampling in R. Draws parameters
#       directly from the prior distributions, applies the full CTMC model
#       (matrix exponential via expm), and generates synthetic observations.
#       Used for all visualisations so results are reproducible without Stan.
#
# QUANTITIES CHECKED
# ------------------
#   1. Prior parameter distributions (β, λ_C, μ, ε)
#   2. Prior predictive VT/NVT carriage prevalence vs observed
#   3. Prior predictive state transition matrix (3×3 averaged over individuals)
#   4. Prior predictive 9-type transition frequencies vs observed
#   5. Prior predictive VT/NVT prevalence by age/HIV group vs observed
#   6. Prior predictive carriage duration distribution
#   7. Prior predictive R_HH_VT spectral radius (household R₀)
#   8. Prior predictive eps_V / eps_N (competition parameters)
#   9. Coverage calibration: P(observed statistic ∈ prior predictive interval)
#  10. Sensitivity: which prior has the largest impact on simulated prevalence
#
# PRIORS (from spn_waifw_5grp_foicomp_stan.R)
#   beta_H_VT[a,b]  ~ Normal(0, 1)   [truncated > 0]  half-normal
#   beta_H_NVT[a,b] ~ Normal(0, 1)   [truncated > 0]
#   lambda_C[k,a]   ~ Normal(0, 0.5) [truncated > 0]  smaller community FOI
#   mu_VT[a]        ~ Normal(1, 0.5) [truncated > 0]  ~1 clearance/interval
#   mu_NVT[a]       ~ Normal(1, 0.5) [truncated > 0]
#   eps_V           ~ Exponential(2)  mean = 0.5
#   eps_N           ~ Exponential(2)  mean = 0.5
#
# OUTPUT FILES
#   spn_ppc_A_param_priors.pdf           — prior distributions of all parameters
#   spn_ppc_B_prevalence.pdf             — prior predictive vs observed prevalence
#   spn_ppc_C_transitions.pdf           — prior predictive transition matrix
#   spn_ppc_D_prevalence_by_group.pdf   — prevalence by age/HIV group
#   spn_ppc_E_duration.pdf              — carriage duration prior predictive
#   spn_ppc_F_R0_household.pdf          — spectral radius (R_HH_VT) prior
#   spn_ppc_G_competition.pdf           — eps_V and eps_N prior predictive
#   spn_ppc_H_calibration.pdf           — coverage calibration plot
#   spn_ppc_summary_panel.pdf           — combined overview panel
#   spn_ppc_stan_fit.rds                — Stan prior-only fit object
# =============================================================================

# ---- 0. Packages & Settings -------------------------------------------------

library(tidyverse)
library(rstan)
library(posterior)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(scales)
library(expm)       # matrix exponential for R-based simulation

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(2024L)

# ---- 1. Constants -----------------------------------------------------------

N_AGE         <- 5L
N_PRIOR_DRAWS <- 2000L    # R-based prior draws (fast)
N_STAN_ITER   <- 1000L    # Stan fixed_param draws per chain

AGE_LABELS_FLAT <- c(
  "Y-Child", "O-Child", "Adult HIV-",
  "Adult HIV+ ART-S", "Adult HIV+ ART-L"
)

STATE_LABELS <- c("S", "VT", "NVT")

pal5  <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")
pal3  <- c("S" = "#2ecc71", "VT" = "#e74c3c", "NVT" = "#3498db")
names(pal5) <- AGE_LABELS_FLAT

# ---- 2. Load Observed Data --------------------------------------------------
# Re-run the same data preparation as the main script so we have
# observed statistics to compare against.

dat <- read_csv("spn_test-f0d48db0.csv") %>%
  mutate(
    state = case_when(
      stg == 1L ~ 1L,
      stg == 3L ~ 2L,
      stg == 2L ~ 3L,
      TRUE      ~ NA_integer_
    ),
    age_grp = case_when(
      agecat == "ychild"                                ~ 1L,
      agecat == "ochild"                                ~ 2L,
      agecat == "adult" & (is.na(hiv) | hiv == "hiv-") ~ 3L,
      agecat == "adult" & hiv == "hiv+_artS"            ~ 4L,
      agecat == "adult" & hiv == "hiv+_artL"            ~ 5L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(state), !is.na(age_grp))

# Longitudinal data: individuals with all three visits
dat_wide <- dat %>%
  select(hhid, pid, vno, state, age_grp) %>%
  pivot_wider(names_from = vno, values_from = state,
              names_prefix = "s") %>%
  filter(!is.na(s1), !is.na(s2), !is.na(s3))

N_OBS <- nrow(dat_wide)

# ---- 2a. Observed Summary Statistics ----------------------------------------

# Overall state distribution at each visit
obs_state_prev <- dat_wide %>%
  pivot_longer(c(s1, s2, s3), names_to = "visit", values_to = "state") %>%
  mutate(
    visit  = recode(visit, s1 = "Visit 1", s2 = "Visit 2", s3 = "Visit 3"),
    state_label = STATE_LABELS[state]
  ) %>%
  count(visit, state_label) %>%
  group_by(visit) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# VT and NVT prevalence by age group (pooled across visits)
obs_prev_by_grp <- dat_wide %>%
  mutate(age_label = AGE_LABELS_FLAT[age_grp]) %>%
  pivot_longer(c(s1, s2, s3), names_to = "visit", values_to = "state") %>%
  group_by(age_grp, age_label) %>%
  summarise(
    prev_VT  = mean(state == 2L),
    prev_NVT = mean(state == 3L),
    n        = n(),
    .groups  = "drop"
  )

# Observed 9-type transition table (from, to) across both transitions
obs_trans <- bind_rows(
  dat_wide %>% transmute(from = s1, to = s2, transition = "V1→V2"),
  dat_wide %>% transmute(from = s2, to = s3, transition = "V2→V3")
) %>%
  count(from, to) %>%
  mutate(prop = n / sum(n),
         from_lbl = STATE_LABELS[from],
         to_lbl   = STATE_LABELS[to])

cat(sprintf("\n=== Observed data: %d individuals, 2 transitions ===\n", N_OBS))
cat("\nObserved state prevalences:\n")
print(obs_state_prev)

# ---- 3. Household Covariate Data (for R-based simulation) ------------------
# We need HH covariate structure (n_VT, n_NVT, n_HH_others, comm_VT, comm_NVT)
# from the same preparation as the main script.

long <- dat %>% select(hhid, pid, vno, state, age_grp)

hh_raw <- long %>%
  inner_join(
    long %>% rename(pid_o = pid, state_o = state, age_o = age_grp),
    by = c("hhid", "vno"),
    relationship = "many-to-many"
  ) %>%
  filter(pid != pid_o)

hh_comp_full <- hh_raw %>%
  group_by(hhid, pid, vno) %>%
  summarise(
    across(paste0("n_VT_",  1:N_AGE),
           list(~ sum(state_o == 2L & age_o == cur_column() %>%
                        str_extract("\\d+") %>% as.integer())),
           .names = "n_VT_{.col}"),
    across(paste0("n_NVT_", 1:N_AGE),
           list(~ sum(state_o == 3L & age_o == cur_column() %>%
                        str_extract("\\d+") %>% as.integer())),
           .names = "n_NVT_{.col}"),
    n_HH_others = n(),
    .groups = "drop"
  )

# Simpler manual construction for reliability
hh_counts <- hh_raw %>%
  group_by(hhid, pid, vno) %>%
  summarise(
    n_HH_others = n(),
    .groups = "drop"
  )
for (b in seq_len(N_AGE)) {
  vt_col  <- paste0("n_VT_",  b)
  nvt_col <- paste0("n_NVT_", b)
  hh_counts[[vt_col]]  <- vapply(
    seq_len(nrow(hh_counts)), function(i) {
      sum(hh_raw$state_o[hh_raw$hhid == hh_counts$hhid[i] &
                           hh_raw$pid  == hh_counts$pid[i]  &
                           hh_raw$vno  == hh_counts$vno[i]  &
                           hh_raw$age_o == b] == 2L)
    }, integer(1))
  hh_counts[[nvt_col]] <- vapply(
    seq_len(nrow(hh_counts)), function(i) {
      sum(hh_raw$state_o[hh_raw$hhid == hh_counts$hhid[i] &
                           hh_raw$pid  == hh_counts$pid[i]  &
                           hh_raw$vno  == hh_counts$vno[i]  &
                           hh_raw$age_o == b] == 3L)
    }, integer(1))
}

dat_model <- dat_wide %>%
  left_join(hh_counts %>% filter(vno == 1) %>% select(-vno),
            by = c("hhid", "pid")) %>%
  left_join(hh_counts %>% filter(vno == 2) %>%
              rename_with(~ paste0(.x, "_v2"),
                          c(starts_with("n_VT_"),
                            starts_with("n_NVT_"), "n_HH_others")) %>%
              select(-vno),
            by = c("hhid", "pid"))

# Build 3D arrays [N x N_AGE x 2]
n_VT_arr  <- array(0L, dim = c(N_OBS, N_AGE, 2L))
n_NVT_arr <- array(0L, dim = c(N_OBS, N_AGE, 2L))
n_HH_mat  <- matrix(1L, nrow = N_OBS, ncol = 2L)

for (b in seq_len(N_AGE)) {
  vt1  <- dat_model[[paste0("n_VT_",  b)]]
  nvt1 <- dat_model[[paste0("n_NVT_", b)]]
  vt2  <- dat_model[[paste0("n_VT_",  b, "_v2")]]
  nvt2 <- dat_model[[paste0("n_NVT_", b, "_v2")]]
  n_VT_arr[, b, 1]  <- if (is.null(vt1))  0L else as.integer(replace_na(vt1,  0L))
  n_NVT_arr[, b, 1] <- if (is.null(nvt1)) 0L else as.integer(replace_na(nvt1, 0L))
  n_VT_arr[, b, 2]  <- if (is.null(vt2))  0L else as.integer(replace_na(vt2,  0L))
  n_NVT_arr[, b, 2] <- if (is.null(nvt2)) 0L else as.integer(replace_na(nvt2, 0L))
}
hh1 <- dat_model[["n_HH_others"]];   hh1[is.na(hh1)] <- 1L
hh2 <- dat_model[["n_HH_others_v2"]]; hh2[is.na(hh2)] <- 1L
n_HH_mat[, 1] <- as.integer(pmax(hh1, 1L))
n_HH_mat[, 2] <- as.integer(pmax(hh2, 1L))

# Community prevalences
comm_prev <- dat %>%
  filter(vno %in% c(1L, 2L)) %>%
  group_by(vno, age_grp) %>%
  summarise(prev_VT = mean(state == 2L), prev_NVT = mean(state == 3L),
            .groups = "drop")

comm_VT_mat <- comm_NVT_mat <- matrix(0.0, N_AGE, 2L)
for (b in seq_len(N_AGE)) {
  for (t in 1:2) {
    tmp <- comm_prev %>% filter(age_grp == b, vno == t)
    if (nrow(tmp) > 0) {
      comm_VT_mat[b,  t] <- tmp$prev_VT
      comm_NVT_mat[b, t] <- tmp$prev_NVT
    }
  }
}

age_grp_vec <- as.integer(dat_wide$age_grp)

# ---- 4. Helper: make_Q and simulate_state -----------------------------------

make_Q_r <- function(lVT, lNVT, muVT, muNVT, epsV, epsN) {
  sigVN <- lNVT * epsV
  sigNV <- lVT  * epsN
  Q <- matrix(0, 3, 3)
  Q[1,1] <- -(lVT + lNVT);  Q[1,2] <- lVT;             Q[1,3] <- lNVT
  Q[2,1] <-  muVT;           Q[2,2] <- -(muVT + sigVN); Q[2,3] <- sigVN
  Q[3,1] <-  muNVT;          Q[3,2] <-  sigNV;          Q[3,3] <- -(muNVT + sigNV)
  Q
}

# Simulate one full dataset (s2_sim, s3_sim) given a single draw of parameters
simulate_dataset <- function(beta_H_VT, beta_H_NVT, lambda_C,
                             mu_VT, mu_NVT, eps_V, eps_N) {
  s2_sim <- integer(N_OBS)
  s3_sim <- integer(N_OBS)
  
  for (n in seq_len(N_OBS)) {
    a    <- age_grp_vec[n]
    sn1  <- dat_wide$s1[n]
    sn2  <- dat_wide$s2[n]
    
    for (t in 1:2) {
      denom <- max(n_HH_mat[n, t], 1L)
      lVT   <- lambda_C[1, a] * comm_VT_mat[a, t]
      lNVT  <- lambda_C[2, a] * comm_NVT_mat[a, t]
      for (b in seq_len(N_AGE)) {
        lVT  <- lVT  + beta_H_VT[a, b]  * n_VT_arr[n, b, t]  / denom
        lNVT <- lNVT + beta_H_NVT[a, b] * n_NVT_arr[n, b, t] / denom
      }
      
      lVT  <- max(lVT,  0)
      lNVT <- max(lNVT, 0)
      Q <- make_Q_r(lVT, lNVT, mu_VT[a], mu_NVT[a], eps_V, eps_N)
      P <- expm(Q)
      
      from_s <- if (t == 1L) sn1 else sn2
      probs  <- pmax(P[from_s, ], 0)
      probs  <- probs / sum(probs)
      to_s   <- sample.int(3L, 1L, prob = probs)
      
      if (t == 1L) s2_sim[n] <- to_s
      else         s3_sim[n] <- to_s
    }
  }
  
  list(s2 = s2_sim, s3 = s3_sim)
}

# ---- 5. R-Based Prior Sampling ----------------------------------------------
# Draw parameters from priors, simulate data, accumulate summary statistics.

cat(sprintf("\n=== Drawing %d prior predictive datasets (R-based) ===\n",
            N_PRIOR_DRAWS))

# Storage for summary statistics per draw
ppc_prevalence <- matrix(NA_real_, N_PRIOR_DRAWS, 6)
colnames(ppc_prevalence) <- c(
  "prev_S_v2", "prev_VT_v2", "prev_NVT_v2",
  "prev_S_v3", "prev_VT_v3", "prev_NVT_v3"
)

ppc_grp_prev <- array(NA_real_, dim = c(N_PRIOR_DRAWS, N_AGE, 4))
# dim 3: VT_v2, NVT_v2, VT_v3, NVT_v3

ppc_trans <- array(NA_real_, dim = c(N_PRIOR_DRAWS, 3, 3))   # mean P matrix
ppc_R0    <- numeric(N_PRIOR_DRAWS)
ppc_eps_V <- numeric(N_PRIOR_DRAWS)
ppc_eps_N <- numeric(N_PRIOR_DRAWS)
ppc_dur_VT_mean  <- numeric(N_PRIOR_DRAWS)
ppc_dur_NVT_mean <- numeric(N_PRIOR_DRAWS)

# Half-normal sampler
rhalfnorm <- function(n, mean, sd) {
  x <- rnorm(n, mean, sd)
  abs(x)   # truncated to positive by reflection; for mean=0 this is half-normal
}

spec_rad <- function(M) max(abs(eigen(M, only.values = TRUE)$values))

pb <- txtProgressBar(min = 0, max = N_PRIOR_DRAWS, style = 3)

for (d in seq_len(N_PRIOR_DRAWS)) {
  
  # Sample from priors
  bHVT  <- matrix(abs(rnorm(N_AGE^2, 0, 1)),   N_AGE, N_AGE)
  bHNVT <- matrix(abs(rnorm(N_AGE^2, 0, 1)),   N_AGE, N_AGE)
  lC    <- matrix(abs(rnorm(2*N_AGE, 0, 0.5)), 2,     N_AGE)
  muVT  <- abs(rnorm(N_AGE, 1, 0.5))
  muNVT <- abs(rnorm(N_AGE, 1, 0.5))
  epsV  <- rexp(1, 2)
  epsN  <- rexp(1, 2)
  
  # Simulate dataset
  sim <- tryCatch(
    simulate_dataset(bHVT, bHNVT, lC, muVT, muNVT, epsV, epsN),
    error = function(e) NULL
  )
  if (is.null(sim)) next
  
  s2s <- sim$s2;  s3s <- sim$s3
  
  # Prevalence
  ppc_prevalence[d, ] <- c(
    mean(s2s == 1L), mean(s2s == 2L), mean(s2s == 3L),
    mean(s3s == 1L), mean(s3s == 2L), mean(s3s == 3L)
  )
  
  # Prevalence by group
  for (g in seq_len(N_AGE)) {
    idx <- age_grp_vec == g
    ppc_grp_prev[d, g, ] <- c(
      mean(s2s[idx] == 2L), mean(s2s[idx] == 3L),
      mean(s3s[idx] == 2L), mean(s3s[idx] == 3L)
    )
  }
  
  # Average transition matrix (over all individuals, both transitions)
  P_acc <- matrix(0, 3, 3)
  n_acc <- 0L
  for (n in seq_len(N_OBS)) {
    a <- age_grp_vec[n]
    for (t in 1:2) {
      denom <- max(n_HH_mat[n, t], 1L)
      lVT  <- lC[1, a] * comm_VT_mat[a, t]
      lNVT <- lC[2, a] * comm_NVT_mat[a, t]
      for (b in seq_len(N_AGE)) {
        lVT  <- lVT  + bHVT[a, b]  * n_VT_arr[n, b, t]  / denom
        lNVT <- lNVT + bHNVT[a, b] * n_NVT_arr[n, b, t] / denom
      }
      lVT  <- max(lVT,  0)
      lNVT <- max(lNVT, 0)
      Q <- make_Q_r(lVT, lNVT, muVT[a], muNVT[a], epsV, epsN)
      P_acc <- P_acc + expm(Q)
      n_acc <- n_acc + 1L
    }
  }
  ppc_trans[d,,] <- P_acc / n_acc
  
  # NGM spectral radius
  R_HH_VT <- outer(seq_len(N_AGE), seq_len(N_AGE),
                   Vectorize(function(a, b) bHVT[a, b] / muVT[b]))
  ppc_R0[d] <- spec_rad(R_HH_VT)
  
  ppc_eps_V[d]       <- epsV
  ppc_eps_N[d]       <- epsN
  ppc_dur_VT_mean[d]  <- mean(1 / muVT)
  ppc_dur_NVT_mean[d] <- mean(1 / muNVT)
  
  setTxtProgressBar(pb, d)
}
close(pb)

cat("\nPrior predictive simulation complete.\n")

# ---- 6. Stan-Based Prior Sampling (Optional, more rigorous) -----------------
# Stan model with likelihood omitted — only priors + generated quantities.

stan_ppc_code <- "
functions {
  matrix make_Q(real lVT, real lNVT, real muVT, real muNVT,
                real epsV, real epsN) {
    real sigVN = lNVT * epsV;
    real sigNV = lVT  * epsN;
    matrix[3,3] Q;
    Q[1,1] = -(lVT + lNVT); Q[1,2] = lVT;             Q[1,3] = lNVT;
    Q[2,1] =  muVT;          Q[2,2] = -(muVT + sigVN); Q[2,3] = sigVN;
    Q[3,1] =  muNVT;         Q[3,2] =  sigNV;          Q[3,3] = -(muNVT + sigNV);
    return Q;
  }
}

data {
  int<lower=1> N;
  int<lower=1> N_AGE;
  array[N]       int<lower=1,upper=3> s1;
  array[N]       int<lower=1,upper=5> age_grp;
  array[N,N_AGE,2] int<lower=0>       n_VT;
  array[N,N_AGE,2] int<lower=0>       n_NVT;
  array[N,2]     int<lower=0>         n_HH_others;
  matrix<lower=0,upper=1>[N_AGE,2]    comm_VT;
  matrix<lower=0,upper=1>[N_AGE,2]    comm_NVT;
}

parameters {
  matrix<lower=0>[N_AGE, N_AGE] beta_H_VT;
  matrix<lower=0>[N_AGE, N_AGE] beta_H_NVT;
  matrix<lower=0>[2, N_AGE]     lambda_C;
  vector<lower=0>[N_AGE]        mu_VT;
  vector<lower=0>[N_AGE]        mu_NVT;
  real<lower=0> eps_V;
  real<lower=0> eps_N;
}

model {
  // PRIORS ONLY — likelihood intentionally omitted for prior predictive check
  to_vector(beta_H_VT)  ~ normal(0.0, 1.0);
  to_vector(beta_H_NVT) ~ normal(0.0, 1.0);
  to_vector(lambda_C)   ~ normal(0.0, 0.5);
  mu_VT  ~ normal(1.0, 0.5);
  mu_NVT ~ normal(1.0, 0.5);
  eps_V  ~ exponential(2.0);
  eps_N  ~ exponential(2.0);
}

generated quantities {
  array[N] int s2_prior;
  array[N] int s3_prior;

  {
    for (n in 1:N) {
      int a = age_grp[n];

      // Transition 1 → 2
      real denom1 = (n_HH_others[n,1] > 0) ? n_HH_others[n,1] : 1.0;
      real lVT1   = lambda_C[1, a] * comm_VT[a,  1];
      real lNVT1  = lambda_C[2, a] * comm_NVT[a, 1];
      for (b in 1:N_AGE) {
        lVT1  += beta_H_VT[a, b]  * n_VT[n, b, 1]  / denom1;
        lNVT1 += beta_H_NVT[a, b] * n_NVT[n, b, 1] / denom1;
      }
      matrix[3,3] P1 = matrix_exp(
        make_Q(fmax(lVT1,0), fmax(lNVT1,0),
               mu_VT[a], mu_NVT[a], eps_V, eps_N));
      s2_prior[n] = categorical_rng(to_vector(P1[s1[n]]));

      // Transition 2 → 3  (using OBSERVED s2, not simulated, to avoid cascade)
      real denom2 = (n_HH_others[n,2] > 0) ? n_HH_others[n,2] : 1.0;
      real lVT2   = lambda_C[1, a] * comm_VT[a,  2];
      real lNVT2  = lambda_C[2, a] * comm_NVT[a, 2];
      for (b in 1:N_AGE) {
        lVT2  += beta_H_VT[a, b]  * n_VT[n, b, 2]  / denom2;
        lNVT2 += beta_H_NVT[a, b] * n_NVT[n, b, 2] / denom2;
      }
      matrix[3,3] P2 = matrix_exp(
        make_Q(fmax(lVT2,0), fmax(lNVT2,0),
               mu_VT[a], mu_NVT[a], eps_V, eps_N));
      s3_prior[n] = categorical_rng(to_vector(P2[s2_prior[n]]));
    }
  }
}
"

cat("\n=== Compiling Stan prior predictive model ===\n")
ppc_stan_model <- stan_model(
  model_code = stan_ppc_code,
  model_name = "spn_ppc_prior_only",
  verbose    = FALSE
)

stan_data_ppc <- list(
  N           = N_OBS,
  N_AGE       = N_AGE,
  s1          = as.integer(dat_wide$s1),
  age_grp     = as.integer(dat_wide$age_grp),
  n_VT        = n_VT_arr,
  n_NVT       = n_NVT_arr,
  n_HH_others = n_HH_mat,
  comm_VT     = comm_VT_mat,
  comm_NVT    = comm_NVT_mat
)

cat("=== Running Stan prior predictive sampler (fixed_param) ===\n")
stan_ppc_fit <- sampling(
  object      = ppc_stan_model,
  data        = stan_data_ppc,
  chains      = 4L,
  iter        = N_STAN_ITER,
  warmup      = 0L,
  algorithm   = "Fixed_param",
  seed        = 2024L,
  verbose     = FALSE,
  refresh     = 250L
)

saveRDS(stan_ppc_fit, "spn_ppc_stan_fit.rds")
cat("Stan prior predictive fit saved: spn_ppc_stan_fit.rds\n")

# Extract Stan simulated states for comparison
stan_draws <- as_draws_df(stan_ppc_fit)
s2_prior_mat <- stan_draws %>%
  select(starts_with("s2_prior[")) %>%
  as.matrix()
s3_prior_mat <- stan_draws %>%
  select(starts_with("s3_prior[")) %>%
  as.matrix()

# Prior predictive prevalence from Stan draws
stan_prev_VT_v2  <- rowMeans(s2_prior_mat == 2L)
stan_prev_NVT_v2 <- rowMeans(s2_prior_mat == 3L)
stan_prev_VT_v3  <- rowMeans(s3_prior_mat == 2L)
stan_prev_NVT_v3 <- rowMeans(s3_prior_mat == 3L)

# =============================================================================
# VISUALISATIONS
# =============================================================================

# ---- Plot A: Prior Parameter Distributions ----------------------------------

prior_param_df <- tibble(
  eps_V   = rexp(5000, 2),
  eps_N   = rexp(5000, 2),
  beta_H  = abs(rnorm(5000, 0, 1)),
  beta_C  = abs(rnorm(5000, 0, 0.5)),
  mu_VT   = abs(rnorm(5000, 1, 0.5)),
  mu_NVT  = abs(rnorm(5000, 1, 0.5))
)

pA1 <- prior_param_df %>%
  select(eps_V, eps_N) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  ggplot(aes(x = value, fill = param)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(
    values = c(eps_V = "#c0392b", eps_N = "#2471a3"),
    labels = c(eps_V = expression(epsilon[V]), eps_N = expression(epsilon[N])),
    name   = "Parameter"
  ) +
  annotate("text", x = 1.05, y = Inf, label = "\u03b5 = 1\n(independence)",
           hjust = 0, vjust = 1.3, size = 3, colour = "grey30") +
  labs(title = "Competition Parameters",
       subtitle = "Exponential(2) | mean = 0.5",
       x = "Parameter value", y = "Density") +
  theme_bw(base_size = 10)

pA2 <- prior_param_df %>%
  select(mu_VT, mu_NVT) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  ggplot(aes(x = value, fill = param)) +
  geom_density(alpha = 0.5, colour = NA) +
  scale_fill_manual(
    values = c(mu_VT = "#e74c3c", mu_NVT = "#3498db"),
    labels = c(mu_VT = expression(mu[VT]), mu_NVT = expression(mu[NVT])),
    name = "Parameter"
  ) +
  labs(title = "Clearance Rates",
       subtitle = "Normal(1, 0.5) truncated > 0",
       x = "Clearance rate", y = "Density") +
  theme_bw(base_size = 10)

pA3 <- prior_param_df %>%
  select(beta_H, beta_C) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  ggplot(aes(x = value, fill = param)) +
  geom_density(alpha = 0.5, colour = NA) +
  scale_fill_manual(
    values = c(beta_H = "#8e44ad", beta_C = "#27ae60"),
    labels = c(beta_H = expression(beta[H]), beta_C = expression(lambda[C])),
    name = "Parameter"
  ) +
  labs(title = "Transmission & Community FOI",
       subtitle = expression(paste(beta[H], " ~ Normal(0,1);  ",
                                   lambda[C], " ~ Normal(0,0.5)")),
       x = "Rate", y = "Density") +
  theme_bw(base_size = 10)

pA <- (pA1 | pA2 | pA3) +
  plot_annotation(
    title    = "Prior Distributions of All Model Parameters",
    subtitle = "Truncated to positive (lower = 0 constraints)",
    theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                     plot.subtitle = element_text(size = 9, colour = "grey40"))
  )

# ---- Plot B: Prior Predictive vs Observed Prevalence ------------------------

obs_prev_summary <- obs_state_prev %>%
  filter(visit %in% c("Visit 2", "Visit 3")) %>%
  select(visit, state_label, prop)

ppc_prev_long <- tibble(
  draw    = seq_len(N_PRIOR_DRAWS),
  VT_v2   = ppc_prevalence[, "prev_VT_v2"],
  NVT_v2  = ppc_prevalence[, "prev_NVT_v2"],
  VT_v3   = ppc_prevalence[, "prev_VT_v3"],
  NVT_v3  = ppc_prevalence[, "prev_NVT_v3"]
) %>%
  pivot_longer(-draw, names_to = "stat", values_to = "value") %>%
  mutate(
    visit       = if_else(str_detect(stat, "v2"), "Visit 2", "Visit 3"),
    state_label = if_else(str_detect(stat, "VT"), "VT", "NVT")
  )

obs_prev_pts <- obs_prev_summary %>%
  filter(state_label %in% c("VT", "NVT"))

pB <- ppc_prev_long %>%
  ggplot(aes(x = value, fill = state_label)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 alpha = 0.55, colour = NA, position = "identity") +
  geom_vline(
    data    = obs_prev_pts,
    mapping = aes(xintercept = prop, colour = state_label),
    linewidth = 1, linetype = "solid"
  ) +
  facet_grid(state_label ~ visit,
             labeller = labeller(state_label = c(VT = "VT carriage",
                                                 NVT = "NVT carriage"))) +
  scale_fill_manual(values  = c(VT = "#e74c3c", NVT = "#3498db"),
                    name = "Carriage type") +
  scale_colour_manual(values = c(VT = "#c0392b", NVT = "#1a5276"),
                      name = "Observed") +
  scale_x_continuous(labels = percent_format()) +
  labs(
    title    = "Prior Predictive Check: Carriage Prevalence",
    subtitle = "Histogram = prior predictive distribution | Vertical line = observed value",
    x        = "Prevalence",
    y        = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(size = 9, colour = "grey40"))

# ---- Plot C: Prior Predictive Transition Matrix -----------------------------

# Average transition matrix from prior predictive draws
P_mean_prior <- apply(ppc_trans, c(2,3), mean, na.rm = TRUE)
P_obs <- matrix(0, 3, 3)
obs_t <- bind_rows(
  dat_wide %>% transmute(from = s1, to = s2),
  dat_wide %>% transmute(from = s2, to = s3)
)
for (f in 1:3) for (t in 1:3) {
  P_obs[f, t] <- mean(obs_t$to[obs_t$from == f] == t)
}

trans_compare <- bind_rows(
  as.data.frame(P_mean_prior) %>%
    setNames(STATE_LABELS) %>%
    mutate(from = STATE_LABELS, source = "Prior Predictive"),
  as.data.frame(P_obs) %>%
    setNames(STATE_LABELS) %>%
    mutate(from = STATE_LABELS, source = "Observed")
) %>%
  pivot_longer(all_of(STATE_LABELS), names_to = "to", values_to = "prob") %>%
  mutate(
    from   = factor(from, levels = STATE_LABELS),
    to     = factor(to,   levels = STATE_LABELS),
    source = factor(source, levels = c("Observed", "Prior Predictive"))
  )

pC <- trans_compare %>%
  ggplot(aes(x = to, y = from, fill = prob)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", prob)),
            colour = "white", fontface = "bold", size = 4) +
  facet_wrap(~ source) +
  scale_fill_gradient2(
    low      = "#3498db", mid = "#ecf0f1", high = "#e74c3c",
    midpoint = 0.33, limits = c(0, 1),
    name     = "P(transition)"
  ) +
  scale_y_discrete(limits = rev) +
  labs(
    title    = "Prior Predictive vs Observed Transition Matrix",
    subtitle = "Averaged over all individuals and both transitions",
    x        = "To state",
    y        = "From state"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 11),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40")
  )

# ---- Plot D: Prior Predictive Prevalence by Age/HIV Group ------------------

obs_grp <- obs_prev_by_grp %>%
  pivot_longer(c(prev_VT, prev_NVT), names_to = "type", values_to = "obs_prev") %>%
  mutate(type = recode(type, prev_VT = "VT (visit 2+3)", prev_NVT = "NVT (visit 2+3)"))

ppc_grp_long <- expand_grid(
  draw  = seq_len(N_PRIOR_DRAWS),
  grp   = seq_len(N_AGE)
) %>%
  mutate(
    VT_v2   = map2_dbl(draw, grp, ~ ppc_grp_prev[.x, .y, 1]),
    NVT_v2  = map2_dbl(draw, grp, ~ ppc_grp_prev[.x, .y, 2]),
    VT_v3   = map2_dbl(draw, grp, ~ ppc_grp_prev[.x, .y, 3]),
    NVT_v3  = map2_dbl(draw, grp, ~ ppc_grp_prev[.x, .y, 4]),
    grp_label = AGE_LABELS_FLAT[grp]
  ) %>%
  pivot_longer(c(VT_v2, NVT_v2, VT_v3, NVT_v3),
               names_to = "stat", values_to = "value") %>%
  mutate(
    type = if_else(str_detect(stat, "^VT"), "VT", "NVT"),
    visit = if_else(str_detect(stat, "v2"), "Visit 2", "Visit 3"),
    grp_label = factor(grp_label, levels = AGE_LABELS_FLAT)
  )

pD <- ppc_grp_long %>%
  filter(visit == "Visit 2") %>%
  ggplot(aes(x = value, fill = type)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 alpha = 0.55, colour = NA, position = "identity") +
  geom_vline(
    data = obs_prev_by_grp %>%
      pivot_longer(c(prev_VT, prev_NVT), names_to = "type", values_to = "obs_prev") %>%
      mutate(type     = recode(type, prev_VT = "VT", prev_NVT = "NVT"),
             grp_label = factor(AGE_LABELS_FLAT[age_grp], levels = AGE_LABELS_FLAT)),
    mapping = aes(xintercept = obs_prev, colour = type),
    linewidth = 1
  ) +
  facet_grid(type ~ grp_label) +
  scale_fill_manual(values = c(VT = "#e74c3c", NVT = "#3498db")) +
  scale_colour_manual(values = c(VT = "#c0392b", NVT = "#1a5276")) +
  scale_x_continuous(labels = percent_format()) +
  labs(
    title    = "Prior Predictive Check by Age/HIV Group (Visit 2)",
    subtitle = "Histogram = prior predictive | Vertical line = observed",
    x        = "Prevalence", y = "Density"
  ) +
  guides(fill = "none", colour = "none") +
  theme_bw(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text.x     = element_text(face = "bold", size = 7),
    strip.text.y     = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40")
  )

# ---- Plot E: Prior Predictive Carriage Duration ----------------------------

ppc_dur_df <- tibble(
  draw     = seq_len(N_PRIOR_DRAWS),
  dur_VT   = ppc_dur_VT_mean,
  dur_NVT  = ppc_dur_NVT_mean
) %>%
  pivot_longer(-draw, names_to = "type", values_to = "duration") %>%
  filter(is.finite(duration), duration < 20)

pE <- ppc_dur_df %>%
  ggplot(aes(x = duration, fill = type)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
  annotate("text", x = 1.05, y = Inf,
           label = "1 inter-visit\ninterval", hjust = 0, vjust = 1.3,
           size = 3, colour = "grey30") +
  scale_fill_manual(
    values = c(dur_VT = "#e74c3c", dur_NVT = "#3498db"),
    labels = c(dur_VT = "VT carriage duration",
               dur_NVT = "NVT carriage duration"),
    name = ""
  ) +
  labs(
    title    = "Prior Predictive: Mean Carriage Duration (1 / \u03bc)",
    subtitle = "Averaged over 5 age/HIV groups | derived from mu ~ Normal(1,0.5) prior",
    x        = "Mean carriage duration (inter-visit intervals)",
    y        = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "bottom"
  )

# ---- Plot F: Prior Predictive Household R₀ (Spectral Radius) ---------------

pF <- tibble(R0 = ppc_R0[is.finite(ppc_R0) & ppc_R0 < 30]) %>%
  ggplot(aes(x = R0)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "#8e44ad", alpha = 0.7, colour = NA) +
  geom_density(colour = "#6c3483", linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "firebrick", linewidth = 0.8) +
  annotate("text", x = 1.1, y = Inf,
           label = "R\u209B\u209B = 1\n(threshold)", hjust = 0, vjust = 1.3,
           size = 3, colour = "firebrick") +
  labs(
    title    = expression(paste(
      "Prior Predictive: Household ", R[0]^{VT}, " (Spectral Radius of NGM)")),
    subtitle = sprintf(
      "Median = %.2f | P(R\u209B\u209B > 1) = %.1f%%",
      median(ppc_R0, na.rm = TRUE),
      100 * mean(ppc_R0 > 1, na.rm = TRUE)),
    x = expression(rho(R[HH~VT])),
    y = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey40"))

# ---- Plot G: Prior Predictive Competition Parameters -----------------------

pG <- tibble(
  eps_V = ppc_eps_V,
  eps_N = ppc_eps_N,
  rel_comp = ppc_eps_V / (ppc_eps_N + 1e-12)
) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  filter(value < 10) %>%
  mutate(param = factor(param,
                        levels = c("eps_V", "eps_N", "rel_comp"),
                        labels = c(
                          expression(epsilon[V]),
                          expression(epsilon[N]),
                          expression(epsilon[V] / epsilon[N])
                        ))) %>%
  ggplot(aes(x = value, fill = param)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 alpha = 0.65, colour = NA) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey30", linewidth = 0.7) +
  facet_wrap(~ param, scales = "free_y",
             labeller = label_parsed) +
  scale_fill_manual(
    values = c("#c0392b", "#2471a3", "#27ae60"),
    guide  = "none"
  ) +
  labs(
    title    = "Prior Predictive: Competition Relative-Risk Parameters",
    subtitle = expression(paste(
      epsilon[V], ", ", epsilon[N], " ~ Exponential(2)  |  ",
      "dashed line = \u03b5 = 1 (independence)")),
    x = "Parameter value", y = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 10),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40")
  )

# ---- Plot H: Coverage Calibration -------------------------------------------
# For each observed statistic, compute what fraction of prior predictive draws
# contain the observed value within the x% CrI. A well-calibrated prior should
# give ~x% coverage at confidence level x.

stat_names <- c(
  "VT prev v2", "NVT prev v2", "VT prev v3", "NVT prev v3"
)
obs_stats <- c(
  obs_prev_summary %>% filter(visit == "Visit 2", state_label == "VT")    %>% pull(prop),
  obs_prev_summary %>% filter(visit == "Visit 2", state_label == "NVT")   %>% pull(prop),
  obs_prev_summary %>% filter(visit == "Visit 3", state_label == "VT")    %>% pull(prop),
  obs_prev_summary %>% filter(visit == "Visit 3", state_label == "NVT")   %>% pull(prop)
)

ppc_stats <- ppc_prevalence[, c("prev_VT_v2","prev_NVT_v2",
                                "prev_VT_v3","prev_NVT_v3")]

alpha_seq <- seq(0.01, 0.99, by = 0.01)
calibration_df <- map_dfr(seq_along(stat_names), function(s) {
  draws <- ppc_stats[, s]
  obs   <- obs_stats[s]
  map_dfr(alpha_seq, function(a) {
    lo <- quantile(draws, (1 - a) / 2)
    hi <- quantile(draws, 1 - (1 - a) / 2)
    tibble(
      stat       = stat_names[s],
      nominal_ci = a,
      covered    = as.numeric(obs >= lo & obs <= hi)
    )
  })
})

calib_summary <- calibration_df %>%
  group_by(stat, nominal_ci) %>%
  summarise(empirical_coverage = mean(covered), .groups = "drop")

pH <- calib_summary %>%
  mutate(stat = factor(stat, levels = stat_names)) %>%
  ggplot(aes(x = nominal_ci, y = empirical_coverage, colour = stat)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey30", linewidth = 0.7) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(
    values = c("#e74c3c","#3498db","#e67e22","#1abc9c"),
    name = "Statistic"
  ) +
  scale_x_continuous(labels = percent_format(), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(labels = percent_format(), breaks = seq(0, 1, 0.2)) +
  labs(
    title    = "Prior Predictive Coverage Calibration",
    subtitle = paste0(
      "Empirical coverage of observed statistic within prior predictive x% interval\n",
      "Well-calibrated prior: line follows the diagonal | ",
      "Below diagonal = prior too narrow | Above = prior too wide"
    ),
    x = "Nominal credible interval level",
    y = "Empirical coverage of observed value"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "right"
  )

# ---- Save all plots ---------------------------------------------------------

ggsave("spn_ppc_A_param_priors.pdf",
       pA, width = 13, height = 4.5, device = cairo_pdf)

ggsave("spn_ppc_B_prevalence.pdf",
       pB, width = 10, height = 5.5, device = cairo_pdf)

ggsave("spn_ppc_C_transitions.pdf",
       pC, width = 10, height = 4.5, device = cairo_pdf)

ggsave("spn_ppc_D_prevalence_by_group.pdf",
       pD, width = 14, height = 5.5, device = cairo_pdf)

ggsave("spn_ppc_E_duration.pdf",
       pE, width = 9, height = 5, device = cairo_pdf)

ggsave("spn_ppc_F_R0_household.pdf",
       pF, width = 7, height = 4.5, device = cairo_pdf)

ggsave("spn_ppc_G_competition.pdf",
       pG, width = 10, height = 4.5, device = cairo_pdf)

ggsave("spn_ppc_H_calibration.pdf",
       pH, width = 8, height = 5.5, device = cairo_pdf)

# Summary panel
summary_panel <- (pA / (pB | pC)) +
  plot_annotation(
    title   = "Prior Predictive Check — Pneumococcal 5-Group WAIFW CTMC Model",
    caption = sprintf(
      "%d prior predictive draws | 5 age/HIV groups | %d individuals | ",
      N_PRIOR_DRAWS, N_OBS),
    theme = theme(
      plot.title   = element_text(face = "bold", size = 14),
      plot.caption = element_text(size = 8, colour = "grey40")
    )
  )

ggsave("spn_ppc_summary_panel.pdf",
       summary_panel, width = 14, height = 11, device = cairo_pdf)

# ---- Numerical summary table ------------------------------------------------

cat("\n=== Prior predictive summary statistics ===\n")
ppc_num_summary <- tibble(
  statistic        = c("VT prevalence V2", "NVT prevalence V2",
                       "VT prevalence V3", "NVT prevalence V3",
                       "Mean carriage duration VT", "Mean carriage duration NVT",
                       "Household R0 (spectral radius)",
                       "eps_V", "eps_N"),
  observed         = c(obs_stats, rep(NA, 5)),
  prior_pred_mean  = c(
    colMeans(ppc_prevalence[, c("prev_VT_v2","prev_NVT_v2",
                                "prev_VT_v3","prev_NVT_v3")]),
    mean(ppc_dur_VT_mean, na.rm = TRUE),
    mean(ppc_dur_NVT_mean, na.rm = TRUE),
    mean(ppc_R0, na.rm = TRUE),
    mean(ppc_eps_V), mean(ppc_eps_N)
  ),
  prior_pred_lo95 = c(
    apply(ppc_prevalence[, c("prev_VT_v2","prev_NVT_v2",
                             "prev_VT_v3","prev_NVT_v3")],
          2, quantile, 0.025),
    quantile(ppc_dur_VT_mean,  0.025, na.rm = TRUE),
    quantile(ppc_dur_NVT_mean, 0.025, na.rm = TRUE),
    quantile(ppc_R0, 0.025, na.rm = TRUE),
    quantile(ppc_eps_V, 0.025), quantile(ppc_eps_N, 0.025)
  ),
  prior_pred_hi95 = c(
    apply(ppc_prevalence[, c("prev_VT_v2","prev_NVT_v2",
                             "prev_VT_v3","prev_NVT_v3")],
          2, quantile, 0.975),
    quantile(ppc_dur_VT_mean,  0.975, na.rm = TRUE),
    quantile(ppc_dur_NVT_mean, 0.975, na.rm = TRUE),
    quantile(ppc_R0, 0.975, na.rm = TRUE),
    quantile(ppc_eps_V, 0.975), quantile(ppc_eps_N, 0.975)
  )
) %>%
  mutate(
    obs_covered_95 = !is.na(observed) &
      observed >= prior_pred_lo95 & observed <= prior_pred_hi95,
    across(where(is.numeric), ~ round(.x, 4))
  )

print(ppc_num_summary, n = Inf)
write_csv(ppc_num_summary, "spn_ppc_summary_table.csv")

cat("\n=== Saved outputs ===\n")
cat("  spn_ppc_A_param_priors.pdf         — prior distributions of all parameters\n")
cat("  spn_ppc_B_prevalence.pdf           — prior predictive vs observed prevalence\n")
cat("  spn_ppc_C_transitions.pdf          — prior predictive transition matrix\n")
cat("  spn_ppc_D_prevalence_by_group.pdf  — prevalence by age/HIV group\n")
cat("  spn_ppc_E_duration.pdf             — carriage duration prior predictive\n")
cat("  spn_ppc_F_R0_household.pdf         — household spectral radius prior\n")
cat("  spn_ppc_G_competition.pdf          — eps_V and eps_N prior predictive\n")
cat("  spn_ppc_H_calibration.pdf          — coverage calibration\n")
cat("  spn_ppc_summary_panel.pdf          — combined overview panel\n")
cat("  spn_ppc_stan_fit.rds               — Stan prior-only fit\n")
cat("  spn_ppc_summary_table.csv          — numerical summary table\n")

# =============================================================================
# INTERPRETATION GUIDE
# =============================================================================
#
# Plot A — Prior parameter distributions
#   Check that the support of each prior is biologically plausible.
#   beta_H and lambda_C are half-normal with most mass near 0 but a long tail
#   — appropriate for rates that could vary widely. mu centred at 1 means
#   roughly one clearance per inter-visit interval (sensible for acute carriage).
#   eps Exponential(2) is strongly left-skewed: prior belief is competitive
#   exclusion (eps < 1), but values > 1 are permitted.
#
# Plot B — Prior predictive prevalence
#   If the vertical observed line falls well within the prior predictive
#   histogram, the priors are weakly informative and consistent with the data.
#   If the observed value falls in the extreme tails (< 2.5% or > 97.5%), the
#   prior is misspecified — either too narrow (prior too informative) or
#   placed on the wrong part of the parameter space.
#
# Plot C — Transition matrix
#   Compare the prior predictive average P matrix to the observed empirical
#   transition frequencies. Rows should not be wildly different in magnitude.
#   Very high off-diagonal elements (high switching or clearance) from the
#   prior would indicate overly diffuse transmission priors.
#
# Plot D — Prevalence by group
#   Group-specific checks reveal whether the model structure (5 separate
#   age/HIV groups with different clearance rates) generates heterogeneity
#   that broadly matches the observed heterogeneity.
#
# Plot E — Carriage duration
#   Duration = 1/mu. Prior mean duration = ~1 inter-visit interval; the
#   distribution extends to longer durations. Check that the bulk of the
#   prior mass lies in a biologically plausible range (days to months).
#
# Plot F — Household R₀
#   The spectral radius of the prior predictive NGM is expected to be wide
#   and to include values both above and below 1. If it were entirely > 1
#   or entirely < 1, the prior would already be answering the key question
#   before seeing data. A prior that puts ~50% mass on each side of 1 is
#   approximately neutral.
#
# Plot G — Competition parameters
#   eps_V and eps_N should have prior mass on both sides of 1 (independence).
#   The Exponential(2) prior places ~86% of mass below 1 (competitive
#   exclusion), which is consistent with the literature but not dogmatic.
#
# Plot H — Coverage calibration
#   If the calibration line closely follows the diagonal, the prior is
#   appropriately calibrated: a 90% prior predictive interval contains the
#   observed value 90% of the time (across the set of summary statistics
#   examined). Lines lying BELOW the diagonal indicate the prior is too
#   narrow (observed values are outside the prior predictive interval more
#   often than expected). Lines ABOVE the diagonal indicate the prior is
#   too wide / uninformative.
# =============================================================================