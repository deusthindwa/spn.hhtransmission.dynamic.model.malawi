#written by Deus
#10/04/2026
#WAIFW household pneumococcal carriage transmission modelling in Malawi

# #====================================================================
# 
# remove.packages(c("StanHeaders", "rstan"))
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

#load packages
pacman::p_load(char = c("lubridate", "tidyverse", "dplyr","here", "rio", "scales", "boot", "magrittr",  "mvtnorm", "zoo", 
                        "patchwork", "ggplotify", "sf", "expm", "INLA", "data.table", "StanHeaders", "bayesplot",
                        "PropCIs", "reshape2","purrr", "msm", "minqa", "optimx", "ggridges", "timetk", "ggbreak", "ggpubr", 
                        "gridExtra", "doParallel", "igraph", "rgdal", "rstan", "expm", "posterior", "expm"))

#====================================================================

#set seed for entire session globally to ensure reproducibility using a task call
#addTaskCallback(function(...) {set.seed(12345);TRUE})

#turn off the global task call for set seed if needed
#removeTaskCallback(1)

#data preparations for analysis
source(here("script", "1_data_preparations.R"))

#description of the study population
source(here("script", "2_pop_characteristics.R"))

#describes the WAIFW household transmission model in Stan langiage
source(here("script", "3_stan_model.R"))


#(TAKES MANY HOURS TO RUN, USE FITTED MODEL OBJECTS IN 'RESULTS' FOLDER TO EXPLORE POSTERIORS)
#compiles and runs the Stan model 
source(here("script", "4_model_run.R"))


#ANALYSIS OF THE POSTERIOR SAMPLES

#generates different estimates and plots from the posterior distribution
source(here("script", "5_posterior_analysis.R"))

#generates different posteriors for sensitivity on missing HIV status
source(here("script", "6_sensitivity_analysis.R"))


# MODEL OVERVIEW
#
#continuous-time Markov chain (CTMC) with 3 carriage states per individual:
#   State 1 — Susceptible (S)  : no carriage
#   State 2 — VT carrier       : vaccine serotype
#   State 3 — NVT carrier      : non-vaccine serotype
#
#five age/HIV groups (WAIFW dimension):
#   Group 1 — ychild       : younger child (0-5 years-old)
#   Group 2 — ochild       : older child (6-17 years-old)
#   Group 3 — adult_hivneg : adult HIV-negative OR unknown HIV status (18+ years-old)
#   Group 4 — adult_hivS   : adult HIV+, short-term ART (hiv+_artS) (18+ years-old)
#   Group 5 — adult_hivL   : adult HIV+, long-term ART  (hiv+_artL) (18+ years-old)
#
#WAIFW 5×5 household transmission matrices:
#   beta_H_VT[a,b]  — rate at which susceptible of group a acquires VT from a household member of group b  (25 pathways)
#   beta_H_NVT[a,b] — same for NVT  (25 pathways)
#
#forces of infection on individual i of group a at transition t:
#   λ_VT(i,t)  = Σ_b beta_H_VT[a,b]  × [n_VT_b_HH  / n_HH_others] + lambda_C[VT,  a]  × prev_VT_community(a, t)
#   λ_NVT(i,t) = Σ_b beta_H_NVT[a,b] × [n_NVT_b_HH / n_HH_others] + lambda_C[NVT, a]  × prev_NVT_community(a, t)
#
#group-specific clearance rates:
#   mu_VT[a]   : VT  clearance rate for carriers of group a
#   mu_NVT[a]  : NVT clearance rate for carriers of group a
#
#FOI-dependent competition
# competition rates (direct serotype switching) are proportional to the force of infection of the competing serotype, modulated by a relative-risk parameter:
#   sigma_VN(i,t) = lambda_NVT(i,t) × eps_V
#   sigma_NV(i,t) = lambda_VT(i,t)  × eps_N
#
#   eps_V : relative risk of acquiring NVT carriage while already carrying VT, compared with a fully susceptible individual.
#           a value < 1 implies partial cross-protection (competitive exclusion)
#           a value = 1 implies independence
#           a value > 1 implies facilitation
#   eps_N : same, for acquiring VT while already carrying NVT
#
#biological interpretation:
#   • sigma_VN(i,t) = 0  whenever lambda_NVT(i,t) = 0: a VT carrier cannot switch to NVT if there is no NVT exposure pressure — the switch requires actual encounter with NVT, not just an intrinsic spontaneous event
#   • Switching rates therefore vary across individuals and over time, proportional to their actual NVT/VT exposure at each transition
#   • eps_V and eps_N are the mechanistically interpretable competition parameters; sigma_VN and sigma_NV are individual-visit-specific derived quantities (not estimated directly).
#
# 3×3 rate matrix Q for individual n of group a at transition t:
#   S   → VT  : lambda_VT(i,t)
#   S   → NVT : lambda_NVT(i,t)
#   VT  → S   : mu_VT[a]
#   VT  → NVT : lambda_NVT(i,t) × eps_V        ←— FOI-dependent
#   NVT → S   : mu_NVT[a]
#   NVT → VT  : lambda_VT(i,t)  × eps_N        ←— FOI-dependent
#
#transition probability matrix: P = expm(Q) over Δt = 1 week
#
#household reproduction number (infector group b's clearance governs duration)
#   R_HH_VT[a,b]  = beta_H_VT[a,b]  / mu_VT[b]
#   R_HH_NVT[a,b] = beta_H_NVT[a,b] / mu_NVT[b]
#
#relative competition:
#   relative_competition = eps_V / eps_N
#   > 1: NVT has relative competitive advantage over VT (eps_V > eps_N)
#   < 1: VT has relative competitive advantage over NVT
#   = 1: symmetric relative susceptibility
#
#total parameters: 72
#   50  β_H  (WAIFW: 5×5 per VT/NVT)
#   10  λ_C  (community FOI: 2 types × 5 groups)
#   10  μ    (clearance:     2 types × 5 groups)
#    2  ε    (eps_V, eps_N)
#
#inference: Hamiltonian Monte Carlo (HMC) in Stan via RStan
#
#data: spn_test-f0d48db0.csv
#   hhid   = household ID        pid    = individual ID
#   vno    = visit number (1,2,3) stg   = carriage state (1=S, 2=VT, 3=NVT)
#   agecat = age category (ychild / ochild / adult)
#   hiv    = HIV/ART status (adult hiv- / hiv+_artS / hiv+_artL)
#   hhsize = household size
#
#summary of all parameters
# beta_H_VT[a,b]      : HH VT  transmission rate; susceptible a from infector
# beta_H_NVT[a,b]     : HH NVT transmission rate; same convention
# lambda_C[1,a]       : community VT  FOI for group
# lambda_C[2,a]       : community NVT FOI for group
# mu_VT[a]            : VT  clearance rate for group
# mu_NVT[a]           : NVT clearance rate for group
# eps_V               : relative risk of NVT acquisition while carrying VT
# eps_N               : relative risk of VT  acquisition while carrying NVT
# [DERIVED] sigma_VN(i,t) = lambda_NVT(i,t) x eps_V  [individual-specific]
# [DERIVED] sigma_NV(i,t) = lambda_VT(i,t)  x eps_N  [individual-specific]
# mean_sigma_VN       : mean(lNVT) x eps_V  [dataset-level mean]
# mean_sigma_NV       : mean(lVT)  x eps_N  [dataset-level mean]
# relative_competition: eps_V / eps_N  (>1 = NVT advantage)
# R_HH_VT[a,b]        : beta_H_VT[a,b]  / mu_VT[b]
# R_HH_NVT[a,b]       : beta_H_NVT[a,b] / mu_NVT[b]
# dur_VT[a]           : 1 / mu_VT[a]  (mean carriage duration)
# dur_NVT[a]          : 1 / mu_NVT[a]
