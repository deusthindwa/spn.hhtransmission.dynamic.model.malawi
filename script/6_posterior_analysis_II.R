#written by Deus
#10/04/2026
#WAIFW household pneumococcal carriage transmission modelling in Malawi

# #====================================================================
#
#decompose the total force of infection (FOI) on each age/HIV group into:
#
#   (A) HOUSEHOLD component: transmission from observed household members
#       lambda_HH_VT(i,t)  = Σ_b  beta_H_VT[a,b]  × n_VT_b(i,t)  / n_HH_others
#       lambda_HH_NVT(i,t) = Σ_b  beta_H_NVT[a,b] × n_NVT_b(i,t) / n_HH_others
#
#   (B) COMMUNITY component: background transmission from outside the household
#       lambda_C_VT(i,t)   = lambda_C[1,a] × pi_VT[a,t]
#       lambda_C_NVT(i,t)  = lambda_C[2,a] × pi_NVT[a,t]
#
#   Total FOI = HH + community  (no interaction between components)

#
#key quantities propagated through the full posterior
# -----------------------------------------------------------
#   1. group-mean HH and community FOI (VT, NVT) per age/HIV group
#   2. HH fraction = lambda_HH / (lambda_HH + lambda_C)
#   3. posterior distributions of HH fraction by group and serotype
#   5. Counterfactual: FOI if household size were fixed (size-adjusted comparison)
#   6. Absolute vs relative attribution across strata
#
# =============================================================================

options(mc.cores = parallel::detectCores())
set.seed(1988L)

N_AGE <- 5L
AGE_LABELS_FLAT <- c("Y-Child", "O-Child", "Adult\nHIV\u2212", "Adult\nHIV+\nART-S", "Adult\nHIV+\nART-L")
AGE_SHORT <- c("YC", "OC", "HIV\u2212", "HIV+S", "HIV+L")
SERO      <- c("VT", "NVT")

pal5 <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")
pal2 <- c(Household = "#e74c3c", Community = "#3498db")
names(pal5) <- AGE_LABELS_FLAT

#load stan posterior draws
fit <- rio::import(here::here("results", "fit_hhbayes_main.rds"))
draws_df <- as_draws_df(fit)
N_DRAWS  <- nrow(draws_df)

#loaded posterior draws
get_matrix_draws <- function(prefix, nrow, ncol) {
  # Returns array [N_DRAWS, nrow, ncol]
  arr <- array(NA_real_, dim = c(N_DRAWS, nrow, ncol))
  for (r in seq_len(nrow)) for (c in seq_len(ncol)) {
    col_name <- paste0(prefix, "[", r, ",", c, "]")
    if (col_name %in% names(draws_df)) {
      arr[, r, c] <- draws_df[[col_name]]
    }
  }
  arr
}

beta_VT_draws  <- get_matrix_draws("beta_H_VT",  N_AGE, N_AGE)
beta_NVT_draws <- get_matrix_draws("beta_H_NVT", N_AGE, N_AGE)
lC_draws       <- get_matrix_draws("lambda_C",   2L,    N_AGE)  # [2 x N_AGE]

mu_VT_draws  <- do.call(cbind, lapply(seq_len(N_AGE), function(a)
  draws_df[[paste0("mu_VT[",  a, "]")]]))
mu_NVT_draws <- do.call(cbind, lapply(seq_len(N_AGE), function(a)
  draws_df[[paste0("mu_NVT[", a, "]")]]))

#-------------------------------------------------------------------------------- 

#reload observed data & covariates
dat <- 
  spn_hhF %>% 
  dplyr::rename('agex' = 'agecat' ) %>%
  dplyr::rename('agecat' = 'agecatx' ) %>%
  dplyr::rename('agecatx' = 'agex' ) %>%
  dplyr::mutate(
    state = case_when(
      stg == 1L ~ 1L,
      stg == 2L ~ 2L,
      stg == 3L ~ 3L,
      TRUE      ~ NA_integer_)) %>%
  dplyr::group_by(hhid) %>%
  dplyr::mutate(hh_hiv_adult = hiv[which(!is.na(hiv) & agecat %in% c("adult"))][1], #extract a single non-missing  HIV value per household (if any)
                hiv = if_else(is.na(hiv) & agecat %in% c("adult"), hh_hiv_adult, hiv)) %>% #fill missing HIV for adults only
  dplyr::ungroup() %>%
  dplyr::select(-hh_hiv_adult) %>%
  dplyr::mutate(
    age_grp = case_when(
      agecat == "ychild"                     ~ 1L,
      agecat == "ochild"                     ~ 2L,
      agecat == "adult" & hiv == "hiv-"      ~ 3L,
      agecat == "adult" & hiv == "hiv+_artS" ~ 4L,
      agecat == "adult" & hiv == "hiv+_artL" ~ 5L,
      TRUE ~ NA_integer_)) %>%
  dplyr::filter(!is.na(state), !is.na(age_grp))

dat_wide <- dat %>%
  select(hhid, pid, vno, state, age_grp) %>%
  pivot_wider(names_from = vno, values_from = state, names_prefix = "s") %>%
  filter(!is.na(s1), !is.na(s2), !is.na(s3))

N_OBS <- nrow(dat_wide)

#-------------------------------------------------------------------------------- 

#community prevalence matrix [N_AGE × 2]
comm_prev <- dat %>%
  filter(vno %in% c(1L, 2L)) %>%
  group_by(vno, age_grp) %>%
  summarise(prev_VT = mean(state == 2L), prev_NVT = mean(state == 3L),
            .groups = "drop")

comm_VT_mat <- comm_NVT_mat <- matrix(0.0, N_AGE, 2L)
for (b in seq_len(N_AGE)) for (t in 1:2) {
  tmp <- comm_prev %>% filter(age_grp == b, vno == t)
  if (nrow(tmp) > 0) {
    comm_VT_mat[b, t]  <- tmp$prev_VT
    comm_NVT_mat[b, t] <- tmp$prev_NVT
  }
}

#household counts [N_OBS × N_AGE × 2] — rebuild from HH raw
long <- dat %>% select(hhid, pid, vno, state, age_grp)
hh_raw <- long %>%
  inner_join(
    long %>% rename(pid_o = pid, state_o = state, age_o = age_grp),
    by = c("hhid", "vno"), relationship = "many-to-many") %>%
  filter(pid != pid_o)

n_VT_arr  <- array(0L, dim = c(N_OBS, N_AGE, 2L))
n_NVT_arr <- array(0L, dim = c(N_OBS, N_AGE, 2L))
n_HH_mat  <- matrix(1L, nrow = N_OBS, ncol = 2L)

for (idx in seq_len(N_OBS)) {
  pid_i  <- dat_wide$pid[idx]
  hhid_i <- dat_wide$hhid[idx]
  for (t in 1:2) {
    sub <- hh_raw %>% filter(hhid == hhid_i, pid == pid_i, vno == t)
    n_HH_mat[idx, t] <- max(nrow(sub), 1L)
    for (b in seq_len(N_AGE)) {
      n_VT_arr[idx, b, t]  <- sum(sub$state_o == 2L & sub$age_o == b)
      n_NVT_arr[idx, b, t] <- sum(sub$state_o == 3L & sub$age_o == b)
    }
  }
}

age_grp_vec <- as.integer(dat_wide$age_grp)

#-------------------------------------------------------------------------------- 

#compute FOI decomposition per posterior draw for 16000 draws
#
# For each draw d, individual n (group a), transition t, compute:
#   lHH_VT[d, n, t]  = Σ_b beta_H_VT[d,a,b] × n_VT[n,b,t] / n_HH[n,t]
#   lC_VT[d, n, t]   = lC[d,1,a] × pi_VT[a,t]
#   lHH_NVT[d, n, t] = similar
#   lC_NVT[d, n, t]  = similar
#
# Then average within group a across individuals and transitions to get
# group-level expected FOI components.
#
# To manage memory with N_DRAWS × N_OBS × 2, we process draws in batches
# and accumulate group-level summaries, not the full individual-level array.
#
#stored per draw × group:
#   lHH_VT_grp[d, a], lC_VT_grp[d, a]
#   lHH_NVT_grp[d, a], lC_NVT_grp[d, a]
#   frac_HH_VT_grp[d, a]  = lHH_VT / (lHH_VT + lC_VT)
#   frac_HH_NVT_grp[d, a] = similar

#0utput arrays
lHH_VT_grp   <- matrix(NA_real_, N_DRAWS, N_AGE)
lC_VT_grp    <- matrix(NA_real_, N_DRAWS, N_AGE)
lHH_NVT_grp  <- matrix(NA_real_, N_DRAWS, N_AGE)
lC_NVT_grp   <- matrix(NA_real_, N_DRAWS, N_AGE)

#WAIFW attribution [d, a_sus, b_inf] — only VT for main analysis
lHH_VT_from_b  <- array(NA_real_, dim = c(N_DRAWS, N_AGE, N_AGE))
lHH_NVT_from_b <- array(NA_real_, dim = c(N_DRAWS, N_AGE, N_AGE))

#per-individual (posterior median) for scatter plot
lHH_VT_ind  <- matrix(NA_real_, N_OBS, 2L)   # [n, t] — posterior median
lC_VT_ind   <- matrix(NA_real_, N_OBS, 2L)
lHH_NVT_ind <- matrix(NA_real_, N_OBS, 2L)
lC_NVT_ind  <- matrix(NA_real_, N_OBS, 2L)

BATCH_SIZE <- 200L
pb <- txtProgressBar(min = 0, max = N_DRAWS, style = 3)

for (d in seq_len(N_DRAWS)) {

  bVT  <- beta_VT_draws[d,,]   # [N_AGE, N_AGE]
  bNVT <- beta_NVT_draws[d,,]
  lC   <- lC_draws[d,,]        # [2, N_AGE]

  grp_lHH_VT  <- numeric(N_AGE)
  grp_lC_VT   <- numeric(N_AGE)
  grp_lHH_NVT <- numeric(N_AGE)
  grp_lC_NVT  <- numeric(N_AGE)
  grp_cnt     <- integer(N_AGE)

  from_b_VT   <- matrix(0, N_AGE, N_AGE)   # [a, b]
  from_b_NVT  <- matrix(0, N_AGE, N_AGE)
  from_b_cnt  <- integer(N_AGE)

  for (n in seq_len(N_OBS)) {
    a <- age_grp_vec[n]
    for (t in 1:2) {
      denom <- max(n_HH_mat[n, t], 1L)

      #community contribution
      lc_vt  <- lC[1, a] * comm_VT_mat[a,  t]
      lc_nvt <- lC[2, a] * comm_NVT_mat[a, t]

      #household contribution — total and by infector group b
      lhh_vt  <- 0
      lhh_nvt <- 0
      for (b in seq_len(N_AGE)) {
        contrib_vt  <- bVT[a, b]  * n_VT_arr[n, b, t]  / denom
        contrib_nvt <- bNVT[a, b] * n_NVT_arr[n, b, t] / denom
        lhh_vt  <- lhh_vt  + contrib_vt
        lhh_nvt <- lhh_nvt + contrib_nvt
        from_b_VT[a, b]  <- from_b_VT[a, b]  + contrib_vt
        from_b_NVT[a, b] <- from_b_NVT[a, b] + contrib_nvt
      }

      grp_lHH_VT[a]  <- grp_lHH_VT[a]  + lhh_vt
      grp_lC_VT[a]   <- grp_lC_VT[a]   + lc_vt
      grp_lHH_NVT[a] <- grp_lHH_NVT[a] + lhh_nvt
      grp_lC_NVT[a]  <- grp_lC_NVT[a]  + lc_nvt
      grp_cnt[a]     <- grp_cnt[a]      + 1L
      from_b_cnt[a]  <- from_b_cnt[a]   + 1L

      #accumulate individual-level (first draw only for scatter plot)
      if (d == 1L) {
        lHH_VT_ind[n, t]  <- lhh_vt
        lC_VT_ind[n, t]   <- lc_vt
        lHH_NVT_ind[n, t] <- lhh_nvt
        lC_NVT_ind[n, t]  <- lc_nvt
      }
    }
  }

  #store group means
  for (a in seq_len(N_AGE)) {
    cnt <- max(grp_cnt[a], 1L)
    lHH_VT_grp[d, a]  <- grp_lHH_VT[a]  / cnt
    lC_VT_grp[d, a]   <- grp_lC_VT[a]   / cnt
    lHH_NVT_grp[d, a] <- grp_lHH_NVT[a] / cnt
    lC_NVT_grp[d, a]  <- grp_lC_NVT[a]  / cnt
    fc <- max(from_b_cnt[a], 1L)
    for (b in seq_len(N_AGE)) {
      lHH_VT_from_b[d, a, b]  <- from_b_VT[a, b]  / fc
      lHH_NVT_from_b[d, a, b] <- from_b_NVT[a, b] / fc
    }
  }

  setTxtProgressBar(pb, d)
}
close(pb)

#-------------------------------------------------------------------------------- 

#compute derived quantities

#HH fraction of total FOI
frac_HH_VT  <- lHH_VT_grp  / (lHH_VT_grp  + lC_VT_grp  + 1e-12)
frac_HH_NVT <- lHH_NVT_grp / (lHH_NVT_grp + lC_NVT_grp + 1e-12)

#posterior median + 50%/95% CrI per group
summarise_grp <- function(mat, label) {
  tibble(
    grp       = seq_len(N_AGE),
    grp_label = AGE_LABELS_FLAT,
    sero      = label,
    mean      = colMeans(mat, na.rm = TRUE),
    median    = apply(mat, 2, median, na.rm = TRUE),
    lo95      = apply(mat, 2, quantile, 0.025, na.rm = TRUE),
    hi95      = apply(mat, 2, quantile, 0.975, na.rm = TRUE),
    lo50      = apply(mat, 2, quantile, 0.25,  na.rm = TRUE),
    hi50      = apply(mat, 2, quantile, 0.75,  na.rm = TRUE)
  )
}

#total FOI of serotype groups
bind_rows(
  summarise_grp(lHH_VT_grp + lC_VT_grp,   "VT")  %>% mutate(source = "Total FOI"),
  summarise_grp(lHH_NVT_grp + lC_NVT_grp,  "NVT") %>% mutate(source = "Total FOI")) %>%
  mutate(grp_label = factor(grp_label, levels = AGE_LABELS_FLAT),
         sero = factor(sero,  levels = c("VT", "NVT")))

#FOI of serotype groups by source
foi_summary <- bind_rows(
  summarise_grp(lHH_VT_grp,   "VT")  %>% mutate(source = "Household"),
  summarise_grp(lC_VT_grp,    "VT")  %>% mutate(source = "Community"),
  summarise_grp(lHH_NVT_grp,  "NVT") %>% mutate(source = "Household"),
  summarise_grp(lC_NVT_grp,   "NVT") %>% mutate(source = "Community")) %>%
  mutate(grp_label = factor(grp_label, levels = AGE_LABELS_FLAT),
         source = factor(source, levels = c("Household", "Community")),
         sero = factor(sero,  levels = c("VT", "NVT")))

#fraction of household FOI
frac_summary <- bind_rows(
  summarise_grp(frac_HH_VT,  "VT"),
  summarise_grp(frac_HH_NVT, "NVT")) %>%
  mutate(grp_label = factor(grp_label, levels = AGE_LABELS_FLAT),
         sero = factor(sero, levels = c("VT", "NVT")))

#posterior median HH fraction of total FOI
frac_summary %>%
  select(grp_label, sero, median, lo95, hi95) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  print(n = Inf)


#VT WAIFW attribution summary [a_sus × b_inf]
waifw_attr_VT <- apply(lHH_VT_from_b, c(2,3), median, na.rm = TRUE)
rownames(waifw_attr_VT) <- paste0("Sus:", AGE_SHORT)
colnames(waifw_attr_VT) <- paste0("Inf:", AGE_SHORT)

#median household VT FOI attribution (susceptible × infector
print(round(waifw_attr_VT, 5))

#NVT WAIFW attribution summary [a_sus × b_inf]
waifw_attr_NVT <- apply(lHH_NVT_from_b, c(2,3), median, na.rm = TRUE)
rownames(waifw_attr_NVT) <- paste0("Sus:", AGE_SHORT)
colnames(waifw_attr_NVT) <- paste0("Inf:", AGE_SHORT)

#median household NVT FOI attribution (susceptible × infector
print(round(waifw_attr_NVT, 5))


#draws long format for ggdist plots
draws_long_grp <- tibble(
  draw     = rep(seq_len(N_DRAWS), N_AGE),
  grp      = rep(seq_len(N_AGE), each = N_DRAWS),
  grp_label = AGE_LABELS_FLAT[grp],
  lHH_VT   = as.vector(lHH_VT_grp),
  lC_VT    = as.vector(lC_VT_grp),
  lHH_NVT  = as.vector(lHH_NVT_grp),
  lC_NVT   = as.vector(lC_NVT_grp),
  fHH_VT   = as.vector(frac_HH_VT),
  fHH_NVT  = as.vector(frac_HH_NVT),
  tot_VT   = lHH_VT + lC_VT,
  tot_NVT  = lHH_NVT + lC_NVT) %>%
  mutate(grp_label = factor(grp_label, levels = AGE_LABELS_FLAT))

#-------------------------------------------------------------------------------- 

#absolute FOI = HH + community side-by-side
pA <- 
  foi_summary %>%
  ggplot(aes(x = grp_label, y = median, colour = source, fill = source, shape = source)) +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.1, linewidth = 1.5, position = position_dodge(width = 0.45), alpha = 0.8) +
  geom_errorbar(aes(ymin = lo50, ymax = hi50), width = 0,   linewidth = 5, position = position_dodge(width = 0.45), alpha = 0.5) +
  geom_point(size = 3, position = position_dodge(width = 0.45), stroke = 1.5) +
  facet_wrap(~ sero, ncol = 2, labeller = labeller(sero = c(VT = "VT carriage", NVT = "NVT carriage"))) +
  scale_colour_manual(values = pal2, name = "FOI source") +
  scale_fill_manual(values   = pal2, name = "FOI source") +
  scale_shape_manual(values  = c(Household = 1, Community = 4), name = "FOI source") +
  scale_y_continuous(limit = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  labs(title = "A", x = "Age/HIV group", y = "Absolute force of infection per week (mean group FOI)") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92"), strip.text = element_text(face = "bold")) +
  theme(plot.title = element_text(face = "bold", size = 20), plot.subtitle = element_text(size = 10, colour = "grey40"), legend.position  = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12), axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12), legend.position = "right") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#stacked proportional attribution bar chart
prop_df <- 
  foi_summary %>%
  group_by(grp_label, sero) %>%
  mutate(total_median = sum(median)) %>%
  ungroup() %>%
  mutate(prop = median / total_median)

pB <- 
  prop_df %>%
  ggplot(aes(x = grp_label, y = prop, fill = source)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  facet_wrap(~ sero, ncol = 2, labeller = labeller(sero = c(VT = "VT FOI", NVT = "NVT FOI"))) +
  scale_fill_manual(values = pal2, name = "FOI source") +
  scale_y_continuous(labels = percent_format(), breaks = seq(0, 1, 0.25)) +
  labs(title = "B", x = "Age/HIV group", y = "Proportion of total force of infection (FOI)") + #Proportional Attribution: Household vs Community FOI by Age/HIV Group | subtitle = "Posterior median | dashed line = equal attribution (50%/50%)"
  theme_bw(base_size = 14) +
  guides(fill = 'none') +
  theme(strip.background = element_rect(fill = "grey92"), strip.text = element_text(face = "bold")) +
  theme(plot.title = element_text(face = "bold", size = 20), plot.subtitle = element_text(size = 10, colour = "grey40"), legend.position  = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12), axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12), legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#posterior density of HH fraction per group
pC <- 
  draws_long_grp %>%
  select(draw, grp_label, fHH_VT, fHH_NVT) %>%
  pivot_longer(c(fHH_VT, fHH_NVT), names_to = "sero", values_to = "frac") %>%
  mutate(sero = recode(sero, fHH_VT  = " VT carriage", fHH_NVT = "NVT carriage")) %>%
  ggplot(aes(x = frac, y = grp_label, fill = grp_label)) +
  stat_halfeye(aes(fill = grp_label), .width = c(0.50, 0.95), point_interval = "median_qi", normalize = "groups", slab_alpha = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  facet_wrap(~ sero, ncol = 2) +
  scale_fill_manual(values = pal5, guide = "none") +
  scale_x_continuous(labels = percent_format(), limits = c(0, 1)) +
  scale_y_discrete(limits = rev) +
  labs(title = "C", x = "Household force of infection (FOI) fraction", y = "Age/HIV group") + #Posterior Distribution of Household FOI Fraction by Group | subtitle = "Proportion of total FOI attributable to household transmission | dashed = 50%", fraction [lambda_HH / (lambda_HH + lambda_C)] 
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92"), strip.text       = element_text(face = "bold")) +
  theme(plot.title = element_text(face = "bold", size = 20), plot.subtitle = element_text(size = 10, colour = "grey40"), legend.position  = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12), axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12), legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

# Combined summary panel
ggsave(here::here("output", "main_CommFOI_2.png"),
       plot = (pA / (pB | pC)) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 18, height = 11, unit = "in", dpi = 300)

#-------------------------------------------------------------------------------- 

#household-size-adjusted FOI comparison
#does the HH/community balance change if we standardise to the same median household size
median_hhsize <- median(rowSums(n_HH_mat) / 2, na.rm = TRUE)

#build HH FOI at median HH size using posterior median beta
beta_VT_med  <- apply(beta_VT_draws,  c(2,3), median)
beta_NVT_med <- apply(beta_NVT_draws, c(2,3), median)
lC_med <- apply(lC_draws, c(2,3), median)

#for each susceptible group a, compute HH FOI as if every other HH member were the same group (each of the 5 groups) at median size:
#this decomposes out individual HH composition effects. we use the actual distribution of HH composition but rescale by median size.
std_hh_df <- 
  expand_grid(sus_grp = seq_len(N_AGE),
              inf_grp = seq_len(N_AGE)) %>%
  mutate(
    #with n_b infectors of type b (we use the group mean) and median HH size:
    std_lHH_VT  = map2_dbl(sus_grp, inf_grp, function(a, b) {
      #mean n_VT_b across all individuals of group a across transitions
      idx <- which(age_grp_vec == a)
      mean_n_vt <- mean(c(n_VT_arr[idx, b, 1], n_VT_arr[idx, b, 2]))
      beta_VT_med[a, b] * mean_n_vt / median_hhsize
    }),
    std_lHH_NVT = map2_dbl(sus_grp, inf_grp, function(a, b) {
      idx <- which(age_grp_vec == a)
      mean_n_nvt <- mean(c(n_NVT_arr[idx, b, 1], n_NVT_arr[idx, b, 2]))
      beta_NVT_med[a, b] * mean_n_nvt / median_hhsize
    }),
    sus_label = factor(AGE_LABELS_FLAT[sus_grp], levels = AGE_LABELS_FLAT)
  ) %>%
  group_by(sus_grp, sus_label) %>%
  summarise(
    std_lHH_VT_total  = sum(std_lHH_VT),
    std_lHH_NVT_total = sum(std_lHH_NVT),
    .groups = "drop"
  )

#community FOI at median community prevalence (visit 1 and 2 average)
comm_lC_df <- 
  tibble(sus_grp   = seq_len(N_AGE),
         sus_label = factor(AGE_LABELS_FLAT, levels = AGE_LABELS_FLAT),
         lC_VT_adj  = lC_med[1, ] * rowMeans(comm_VT_mat),
         lC_NVT_adj = lC_med[2, ] * rowMeans(comm_NVT_mat))

std_compare <- 
  std_hh_df %>%
  left_join(comm_lC_df, by = c("sus_grp", "sus_label")) %>%
  pivot_longer(cols = c(std_lHH_VT_total, std_lHH_NVT_total, lC_VT_adj, lC_NVT_adj), names_to  = "stat", values_to = "foi") %>%
  mutate(sero   = if_else(str_detect(stat, "NVT"),  "NVT",  "VT"), 
         source = if_else(str_detect(stat, "HH"),  "Household (size-adjusted)", "Community"),
         sero   = factor(sero,   levels = c("VT", "NVT")),
         source = factor(source, levels = c("Community", "Household (size-adjusted)")))

pD <- 
  std_compare %>%
  ggplot(aes(x = sus_label, y = foi, fill = source)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  facet_wrap(~ sero, ncol = 2, labeller = labeller(sero = c(VT = "VT FOI", NVT = "NVT FOI"))) +
  scale_fill_manual(values = c("Household (size-adjusted)" = "#c0392b", "Community" = "#2980b9"), name = "") +
  labs(title = '', x = "Susceptible age/HIV group", y = "Standardised FOI (posterior median \u03b2 \u00D7 mean HH counts / median HH size)") +  
  theme_bw(base_size = 11) +
  scale_y_continuous(limit = c(0, 0.75), breaks = seq(0, 0.75, 0.1)) +
  theme(strip.background = element_rect(fill = "grey92"), strip.text = element_text(face = "bold")) +
  theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9, colour = "grey40"), legend.position  = "bottom") +
  theme(plot.title = element_text(face = "bold", size = 20), plot.subtitle = element_text(size = 10, colour = "grey40"), legend.position  = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12), axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12), legend.position = "bottom") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggsave(here::here("output", "suppl_HHsize.png"),
       plot = pD,
       width = 12, height = 8, unit = "in", dpi = 300)

#-------------------------------------------------------------------------------- 

# INTERPRETATION GUIDE
#
# PLOT A — Absolute FOI
#   Shows raw posterior magnitude of HH and community FOI per group.
#   A group where the HH credible interval lies entirely above the community
#   interval is one where household contact is the dominant route of
#   acquisition. Use this to prioritise household-targeted interventions.
#
# PLOT B — Proportional stacked bars
#   Normalises Plot A so both routes sum to 100% per group.
#   Directly reads as "X% of this group's VT exposure comes from the household".
#   Groups above 50% HH are primarily household-driven; those below are mainly
#   community-driven. Differences across HIV status in adults reveal whether
#   HIV-related immune changes alter the relative importance of each route.
#
# PLOT C — Posterior density of HH fraction
#   Full posterior uncertainty in the HH fraction. If the 95% CrI spans both
#   sides of 50%, the data cannot distinguish whether household or community
#   is the dominant source for that group. Wide CrIs indicate where more data
#   or stronger priors are needed.
#
# PLOT D — Household-size-adjusted
#   Removes compositional effects of different household sizes.
#   Standardises to the median observed HH size, making groups comparable
#   on equal footing. If the HH:community ratio changes substantially from
#   Plot B to Plot C, variation in household size is an important confounder.
#
# POLICY RELEVANCE
#   Groups with high HH fractions benefit most from household-level
#   interventions (contact tracing, household vaccination, decolonisation).
#   Groups where community FOI dominates need population-wide approaches.
#   The HIV+ ART groups are of special interest: if their HH fraction is
#   elevated, it may reflect that close household carers drive their exposure,
#   making household vaccination of carers particularly beneficial.
