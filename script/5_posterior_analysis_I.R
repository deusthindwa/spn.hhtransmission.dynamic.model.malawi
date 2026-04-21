#written by Deus
#10/04/2026
#WAIFW household pneumococcal carriage transmission modelling in Malawi

# #====================================================================

#improt the posterior dataset
fit <- rio::import(here::here("results", "fit_hhbayes_main.rds"))

#convergence Diagnostics
cat("\n=== HMC diagnostics ===\n")
check_hmc_diagnostics(fit)

#extract posterior draws
post_draws <- as_draws_array(fit)
post_df    <- as_draws_df(fit)

#summarise posterior draws
rhat_tbl <- 
  posterior::summarise_draws(post_draws, 
                             mean, sd,
                             ~ quantile(.x, c(0.025, 0.50, 0.975)),rhat, ess_bulk, ess_tail)
rhat_tbl

#list all inferred parameters
betaVT_params  <- c(paste0("beta_H_VT[",  1:5,",",1, "]"), 
                  paste0("beta_H_VT[", 1:5, ",",2, "]"),
                  paste0("beta_H_VT[", 1:5, ",",3, "]"),
                  paste0("beta_H_VT[", 1:5, ",",4, "]"),
                  paste0("beta_H_VT[", 1:5, ",",5, "]"))
                  
betaNVT_params<- c(paste0("beta_H_NVT[",  1:5,",",1, "]"), 
                  paste0("beta_H_NVT[", 1:5, ",",2, "]"),
                  paste0("beta_H_NVT[", 1:5, ",",3, "]"),
                  paste0("beta_H_NVT[", 1:5, ",",4, "]"),
                  paste0("beta_H_NVT[", 1:5, ",",5, "]"))

lambda_params<- c(paste0("lambda_C[", 1:2,", ",1, "]"), 
                  paste0("lambda_C[", 1:2, ",",2, "]"),
                  paste0("lambda_C[", 1:2, ",",3, "]"),
                  paste0("lambda_C[", 1:2, ",",4, "]"),
                  paste0("lambda_C[", 1:2, ",",5, "]"))

mu_params  <- c(paste0("mu_VT[",  1:5, "]"), paste0("mu_NVT[", 1:5, "]"))

eps_params <- c("eps_V", "eps_N", "mean_sigma_VN", "mean_sigma_NV") #

#group-specific within household transmission rates (R-hat + ESS)
print(rhat_tbl %>%
        filter(variable %in% betaVT_params) %>%
        mutate(group_to = rep(AGE_LABELS_FLAT, 5), 
               type  = rep(c("VT"), each = 25)) %>%
        select(type, group_to, mean, sd, `2.5%`, `50%`, `97.5%`, rhat, ess_bulk))

print(rhat_tbl %>%
        filter(variable %in% betaNVT_params) %>%
        mutate(group_to = rep(AGE_LABELS_FLAT, 5), 
               type  = rep(c("NVT"), each = 25)) %>%
        select(type, group_to, mean, sd, `2.5%`, `50%`, `97.5%`, rhat, ess_bulk))

#group-specific community force of infection (R-hat + ESS)
print(rhat_tbl %>%
        filter(variable %in% mu_params) %>%
        mutate(group = rep(AGE_LABELS_FLAT, 2), type  = rep(c("VT","NVT"), each = 5)) %>%
        select(type, group, mean, sd, `2.5%`, `50%`, `97.5%`, rhat, ess_bulk))

#group-specific clearance rates (R-hat + ESS)
print(rhat_tbl %>%
        filter(variable %in% mu_params) %>%
        mutate(group = rep(AGE_LABELS_FLAT, 2), type  = rep(c("VT","NVT"), each = 5)) %>%
        select(type, group, mean, sd, `2.5%`, `50%`, `97.5%`, rhat, ess_bulk))

#competition parameters (R-hat + ESS)
print(rhat_tbl %>%
        filter(variable %in% eps_params) %>%
        select(variable, mean, sd, `2.5%`, `50%`, `97.5%`, rhat, ess_bulk)
)

#------------------------------------------------------------------------------

#trace and pairs plots
#transmissition
ggsave(here::here("output", "suppl_p_densp_betaV.png"), plot = bayesplot::mcmc_pairs(post_draws, pars = betaVT_params, diag_fun = "dens"), width = 30, height = 20, unit = "in", dpi = 300)

ggsave(here::here("output", "suppl_p_densp_betaN.png"), plot = bayesplot::mcmc_pairs(post_draws, pars = betaNVT_params, diag_fun = "dens"), width = 30, height = 20, unit = "in", dpi = 300)

#foi
ggsave(here::here("output", "suppl_p_densp_foi.png"), plot = bayesplot::mcmc_pairs(post_draws, pars = lambda_params, diag_fun = "dens"), width = 15, height = 10, unit = "in", dpi = 300)

#clearance
ggsave(here::here("output", "suppl_p_densp_mu.png"), plot = bayesplot::mcmc_pairs(post_draws, pars = mu_params, diag_fun = "dens"), width = 15, height = 10, unit = "in", dpi = 300)

#competition
ggsave(here::here("output", "suppl_p_densp_eps.png"), plot = bayesplot::mcmc_pairs(post_draws, pars = eps_params, diag_fun = "dens"), width = 15, height = 10, unit = "in", dpi = 300)

#------------------------------------------------------------------------------

#WAIFW 5x5 Heatmap Helper function
summarise_matrix_param <- function(prefix, n_age = N_AGE) {
  expand.grid(a = seq_len(n_age), b = seq_len(n_age), KEEP.OUT.ATTRS = FALSE) %>%
    mutate(param = sprintf("%s[%d,%d]", prefix, a, b),
           susc_label = factor(AGE_LABELS[a], levels = AGE_LABELS),
           inf_label  = factor(AGE_LABELS[b], levels = AGE_LABELS)) %>%
    rowwise() %>%
    mutate(med  = median(post_df[[param]]),
           lo95 = quantile(post_df[[param]], 0.025),
           hi95 = quantile(post_df[[param]], 0.975)) %>%
    ungroup()
}

plot_waifw_5x5 <- function(df, title, fill_label = expression(beta), text_size = 2.8) {
  ggplot(df, aes(x = inf_label, y = susc_label, fill = med)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f\n[%.2f,%.2f]", med,lo95, hi95)), size = text_size, colour = "black",  lineheight = 1.1) +
    scale_fill_gradient2(low = "gray90", mid = "orange", high = "red", midpoint = median(df$med), name = fill_label) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(AGE_LABELS)) +
    labs(title    = title, x = "Infector group (who transmits)", y = "Susceptible group (who acquires)" ) + #subtitle = "Posterior median (95% CrI) | Rows = susceptible, Cols = infector",
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 10, face = "bold"), axis.text.x = element_text(angle = 0, hjust = 0.5), plot.title = element_text(face = "bold", size = 20), plot.subtitle = element_text(size = 10, colour = "grey40"), legend.position = "right") +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
}

#WAIFW Beta Heatmaps
waifw_VT  <- summarise_matrix_param("beta_H_VT")
waifw_NVT <- summarise_matrix_param("beta_H_NVT")
betaHH <- print(plot_waifw_5x5(waifw_VT,  "A") + labs(subtitle = "Household VT Transmission (WAIFW \u03b2)") | 
                plot_waifw_5x5(waifw_NVT, "C") + labs(subtitle = "Household NVT Transmission (WAIFW \u03b2)")) 

#R_HH Heatmaps
R_VT  <- summarise_matrix_param("R_HH_VT")
R_NVT <- summarise_matrix_param("R_HH_NVT")
repHH <- print(plot_waifw_5x5(R_VT,  "B", fill_label = expression(R[HH])) + labs(subtitle = "Household VT reproduction number (R)") | 
               plot_waifw_5x5(R_NVT, "D", fill_label = expression(R[HH])) + labs(subtitle = "Household NVT reproduction number (R)"))

#combined plot
WAIFW_main <- betaHH/repHH

#combined plots
ggsave(here::here("output", "main_WAIFW.png"), 
       plot = WAIFW_main, 
       width = 15, 
       height = 10, 
       unit = "in", 
       dpi = 300)

#------------------------------------------------------------------------------

#community FOI
lambda_C_tbl <- expand.grid(type_idx = 1:2, age_idx = seq_len(N_AGE), KEEP.OUT.ATTRS = FALSE) %>%
  mutate(param     = sprintf("lambda_C[%d,%d]", type_idx, age_idx),
         type      = c("VT","NVT")[type_idx],
         age_label = factor(AGE_LABELS_FLAT[age_idx], levels = AGE_LABELS_FLAT)) %>%
  rowwise() %>%
  mutate(med  = median(post_df[[param]]),
         lo95 = quantile(post_df[[param]], 0.025),
         hi95 = quantile(post_df[[param]], 0.975)) %>%
  ungroup()

p_comm <- 
  ggplot(lambda_C_tbl, aes(x = age_label, y = med, colour = type, group = type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.5) +
  geom_point(size = 3.5, position = position_dodge(0.35), shape = 21, stroke = 2) +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.1, linewidth = 0.9, position = position_dodge(0.35)) +
  scale_colour_manual(values = c(VT = "#c0392b", NVT = "#2980b9"), name = "Serotype") +
  labs(x = "Age/HIV group", y = expression(lambda^'C'~(Community~acquisition~per~week ))) + 
  scale_y_continuous(limit = c(0, 1.3), breaks = seq(0, 1.3, 0.2)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, size = 12), axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12), plot.title = element_text(face = "bold"), legend.position = 'right') +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

print(p_comm)

#combined plots
ggsave(here::here("output", "main_CommFOI.png"),
       plot = p_comm,
       width = 10, height = 5, unit = "in", dpi = 300)

#------------------------------------------------------------------------------

#group-specific clearance rate posteriors
mu_tbl <- 
  bind_rows(
    tibble(param = paste0("mu_VT[",  1:5, "]"),
           type  = "VT",
           group = factor(AGE_LABELS_FLAT, levels = AGE_LABELS_FLAT)),
    tibble(param = paste0("mu_NVT[", 1:5, "]"),
           type  = "NVT",
           group = factor(AGE_LABELS_FLAT, levels = AGE_LABELS_FLAT))) %>%
  rowwise() %>%
  mutate(med  = median(post_df[[param]]),
         lo95 = quantile(post_df[[param]], 0.025),
         hi95 = quantile(post_df[[param]], 0.975),
         lo50 = quantile(post_df[[param]], 0.25),
         hi50 = quantile(post_df[[param]], 0.75)) %>%
  ungroup()

p_mu <- 
  ggplot(mu_tbl, aes(x = group, y = med, colour = type, group = type)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 1.5, position = position_dodge(0.45), alpha = 0.8) +
  geom_linerange(aes(ymin = lo50, ymax = hi50), linewidth = 5, position = position_dodge(0.45), alpha = 0.5) +
  geom_point(size = 3.5, position = position_dodge(0.45)) +
  scale_colour_manual(values = c(VT = "#c0392b", NVT = "#2980b9"), name = "Serotype group") +
  labs(title = "A", x = "Age/HIV", y = expression(mu~(Carriage~clearance~per~week))) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(limit = c(0, 1.6), breaks = seq(0, 1.6, 0.2)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12), axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12), plot.title = element_text(face = "bold", size = 20), legend.position = "right") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

print(p_mu)


dur_tbl <- 
  mu_tbl %>%
  mutate(dur_med  = 1 / med,
         dur_lo95 = 1 / hi95,
         dur_hi95 = 1 / lo95,
         dur_lo50 = 1 / hi50,
         dur_hi50 = 1 / lo50)

p_dur <- 
  ggplot(dur_tbl, aes(x = group, y = dur_med, colour = type, group = type)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_linerange(aes(ymin = dur_lo95, ymax = dur_hi95), linewidth = 1.0, position = position_dodge(0.45), alpha = 0.45) +
  geom_linerange(aes(ymin = dur_lo50, ymax = dur_hi50), linewidth = 2.5, position = position_dodge(0.45), alpha = 0.7) +
  geom_point(size = 3.5, position = position_dodge(0.45)) +
  coord_cartesian(ylim = c(0,20)) +
  scale_colour_manual(values = c(VT = "#c0392b", NVT = "#2980b9"), name = "Serotype") +
  labs(x = "Age / HIV group", y = "Mean duration (Weeks)") + #    # title    = "Mean Carriage Duration (1/\u03bc) by Group",# subtitle = "Longer duration in HIV+ groups reflects immune suppression",
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12), axis.text.y = element_text(angle = 0, hjust = 1, size = 12), plot.title = element_text(face = "bold"), legend.position = "top") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

print(p_dur)


#------------------------------------------------------------------------------

#competition parameter posteriors

#posterior densities of eps_V and eps_N
p_eps_dens <- 
  mcmc_areas(post_draws,
             pars = c("eps_V", "eps_N"),
             prob = 0.80,
             prob_outer = 0.95) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "firebrick", linewidth = 0.8) +
  annotate("text", x = 1.02, y = 2.5, label = "eps = 1\n(independence)", hjust = 0, size = 3, colour = "firebrick") +
  labs(title = "Posterior Densities — Relative-Risk Competition Parameters",
       subtitle = paste0("\u03b5_V: relative risk of NVT acquisition while carrying VT  |  ","\u03b5_N: relative risk of VT acquisition while carrying NVT\n",
                         "Dashed line: eps = 1 (no interaction). eps < 1 = competitive exclusion; ",
                         "eps > 1 = facilitation" )) +
  theme_minimal(base_size = 20) +
  theme(plot.subtitle = element_text(size = 8, colour = "grey40"))

print(p_eps_dens)

#point estimates + 95% CrI for eps_V and eps_N, with reference line
eps_summary <- 
  tibble(param = c("eps_V", "eps_N"),
         label = c(expression(epsilon[V]~"(VT carrier \u2192 NVT risk)"), 
                   expression(epsilon[N]~"(NVT carrier \u2192 VT risk)")),
         label_text = c("εV (VT carrier \u2192 NVT risk)", # \u03b5\ \u1D5B
                        "εN (NVT carrier \u2192 VT risk)")) %>% #u2099
  
  rowwise() %>%
  mutate(med  = median(post_df[[param]]),
         lo95 = quantile(post_df[[param]], 0.025),
         hi95 = quantile(post_df[[param]], 0.975),
         lo50 = quantile(post_df[[param]], 0.25),
         hi50 = quantile(post_df[[param]], 0.75)) %>%
  ungroup()

p_eps_point <- 
  ggplot(eps_summary, aes(x = label_text, y = med, colour = label_text)) +
  #geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick", linewidth = 0.8) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 2.0, alpha = 0.8) +
  geom_linerange(aes(ymin = lo50, ymax = hi50), linewidth = 5.0, alpha = 0.5) +
  geom_point(size = 5) +
  scale_y_continuous(limit = c(0, 0.85), breaks = seq(0, 0.8, 0.2)) +
  scale_colour_manual(values = c("εV (VT carrier \u2192 NVT risk)"  = "#8e44ad", "εN (NVT carrier \u2192 VT risk)" = "#16a085"), guide = "none") +
  annotate("text", x = 0.55, y = 1.03, label = "No interaction", hjust = 0, size = 3.2, colour = "firebrick", fontface = "italic") +
  labs(title = "B", x = "Competition susceptibility parameter", y = expression(epsilon~"(Relative risk)")) + #subtitle = "Thick bar: 50% CrI | Thin bar: 95% CrI | Dashed: eps = 1", Competition Relative-Risk Estimates
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12), axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12), plot.title = element_text(face = "bold", size = 20)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

print(p_eps_point)

#relative competition summary
rc_draws <- post_df$relative_competition

cat(sprintf("\nrelative competition (eps_V / eps_N)\n  Median:   %.3f\n  95%% CrI: [%.3f, %.3f]\n  >> %s\n",
            median(rc_draws),
            quantile(rc_draws, 0.025),
            quantile(rc_draws, 0.975),
            ifelse(median(rc_draws) > 1, "NVT has relative competitive advantage (eps_V > eps_N)", "VT has relative competitive advantage (eps_V < eps_N)")))

#mean realised sigma values (at dataset-level mean FOI)
cat("\nrealised mean competition rates at mean FOI\n")
cat(sprintf("mean_sigma_VN = mean(lNVT) x eps_V:  median = %.4f  [%.4f, %.4f]\n",
            median(post_df$mean_sigma_VN),
            quantile(post_df$mean_sigma_VN, 0.025),
            quantile(post_df$mean_sigma_VN, 0.975)))

cat(sprintf("mean_sigma_NV = mean(lVT)  x eps_N:  median = %.4f  [%.4f, %.4f]\n",
            median(post_df$mean_sigma_NV),
            quantile(post_df$mean_sigma_NV, 0.025),
            quantile(post_df$mean_sigma_NV, 0.975)))

#------------------------------------------------------------------------------

#sensitivity of sigma to FOI — show how sigma_VN varies across a range of plausible lNVT values given the posterior of eps_V
eps_V_draws  <- post_df$eps_V
eps_N_draws  <- post_df$eps_N
lNVT_grid    <- seq(0, 0.5, by = 0.01)
lVT_grid     <- seq(0, 0.5, by = 0.01)

sigma_VN_grid <- outer(lNVT_grid, eps_V_draws, "*")   #[grid x draws]
sigma_NV_grid <- outer(lVT_grid,  eps_N_draws, "*")

sigma_VN_df <- tibble(
  lNVT      = lNVT_grid,
  med       = apply(sigma_VN_grid, 1, median),
  lo95      = apply(sigma_VN_grid, 1, quantile, 0.025),
  hi95      = apply(sigma_VN_grid, 1, quantile, 0.975),
  type      = "ϕVN = λNVT x εV")

sigma_NV_df <- tibble(
  lNVT      = lVT_grid,
  med       = apply(sigma_NV_grid, 1, median),
  lo95      = apply(sigma_NV_grid, 1, quantile, 0.025),
  hi95      = apply(sigma_NV_grid, 1, quantile, 0.975),
  type      = "ϕNV = λVT x εN")

p_sigma_curve <- 
  bind_rows(sigma_VN_df, sigma_NV_df) %>%
  ggplot(aes(x = lNVT, y = med, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  scale_colour_manual(values = c("ϕVN = λNVT x εV" = "#8e44ad", "ϕNV = λVT x εN"  = "#16a085"), name = NULL) +
  scale_fill_manual(values = c("ϕVN = λNVT x εV" = "#8e44ad", "ϕNV = λVT x εN"  = "#16a085"), guide = "none") +
  annotate("text", x = 0.42, y = 0.43, label = "ε=1 (reference)", size = 3, colour = "grey40", fontface = "italic") +
  labs(title = "C", x = "λ (Total FOI of competing serotype per week)", y = "ϕ (Switching rate per week)") + #Realised Competition Rate as a Function of FOI | subtitle = paste0("Shaded: 95% posterior CrI | Dashed: sigma = FOI (eps = 1 reference)\n", "x-axis: force of infection of the competing serotype")
  theme_minimal(base_size = 14) +
  scale_x_continuous(limit = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  scale_y_continuous(limit = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12), axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12), plot.title = element_text(face = "bold", size = 20), legend.position = "right") +
  theme(plot.subtitle = element_text(size = 8.5, colour = "grey40")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

print(p_sigma_curve)


#combined plots
ggsave(here::here("output", "main_ClearComp.png"),
       plot = p_mu / (p_eps_point | p_sigma_curve),
       width = 15, height = 9, unit = "in", dpi = 300)

#------------------------------------------------------------------------------

#posterior predictive checks
s2_rep_mat <- as.matrix(fit, pars = "s2_rep")
s3_rep_mat <- as.matrix(fit, pars = "s3_rep")

compute_ppc_df <- function(rep_mat, obs_vec, visit_label) {
  rep_props <- apply(rep_mat, 1, function(row)
    table(factor(row, 1:3)) / length(row))
  tibble(state = 1:3,
         state_label = c("Susceptible","VT","NVT"),
         obs = as.numeric(table(factor(obs_vec, 1:3))) / length(obs_vec),
         pred_med = apply(rep_props, 1, median),
         pred_lo95 = apply(rep_props, 1, quantile, 0.025),
         pred_hi95 = apply(rep_props, 1, quantile, 0.975),
         visit = visit_label
  )
}

ppc_df <- 
  bind_rows(compute_ppc_df(s2_rep_mat, dat_model$s2, "Visit 2  (v1\u2192v2)"),
            compute_ppc_df(s3_rep_mat, dat_model$s3, "Visit 3  (v2\u2192v3)")) %>%
  dplyr::mutate(state_label = factor(state_label, levels = c('Susceptible', 'NVT', 'VT')))

p_ppc <- 
  ggplot(ppc_df, aes(x = state_label, fill = state_label)) +
  geom_col(aes(y = pred_med), alpha = 0.65, colour = NA) +
  geom_errorbar(aes(ymin = pred_lo95, ymax = pred_hi95), width = 0.1, linewidth = 0.85) +
  geom_point(aes(y = obs), shape = 21, size = 4, fill = "white", colour = "black", stroke = 1.5) +
  facet_wrap(~ visit) +
  scale_fill_manual(values = c(Susceptible = "#95a5a6", VT = "#c0392b", NVT = "#2980b9"), guide  = "none") +
  labs(title = "", x = "Pneumococcal carriage state", y = "Proportion") + #Posterior Predictive Check subtitle = "Bars: predicted median (95% CrI)  |  Points: observed proportions",
  theme_minimal(base_size = 14) +
  scale_y_continuous(limit = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  theme(strip.text  = element_text(face = "bold"), plot.title  = element_text(face = "bold"), axis.text.x = element_text(angle = 0, hjust = 0.5, size =12), axis.text.y = element_text(angle = 0, hjust = 0.5, size =12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

print(p_ppc)

#combined plots
ggsave(here::here("output", "main_PosterioChecks.png"),
       plot = p_ppc,
       width = 10, height = 6, unit = "in", dpi = 300)
