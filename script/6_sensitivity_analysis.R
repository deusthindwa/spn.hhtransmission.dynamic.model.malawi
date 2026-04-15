#written by Deus
#10/04/2026
#WAIFW household pneumococcal carriage transmission modelling in Malawi

#====================================================================

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
    theme(panel.grid = element_blank(), axis.text = element_text(size = 10, face = "bold"), axis.text.x   = element_text(angle = 0, hjust = 0.5)) + 
    theme(plot.title = element_text(face = "bold", size = 20), plot.subtitle = element_text(size = 10, colour = "grey40"), legend.position = "right") +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
}


#====================================================================
  
#import the posterior dataset
fit_main <- rio::import(here::here("results", "fit_hhbayes_main.rds"))

#extract posterior draws
post_draws <- as_draws_array(fit_main)
post_df    <- as_draws_df(fit_main)

#summarise posterior draws
rhat_tbl <- posterior::summarise_draws(post_draws, mean, sd, ~ quantile(.x, c(0.025, 0.50, 0.975)),rhat, ess_bulk, ess_tail)
rhat_tbl

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

#WAIFW Beta Heatmaps
waifw_VT_fit_main  <- summarise_matrix_param("beta_H_VT")
waifw_NVT_fit_main <- summarise_matrix_param("beta_H_NVT")
  
#====================================================================

#import the posterior dataset
fit_artL <- rio::import(here::here("results", "fit_hhbayes_longART.rds"))

#extract posterior draws
post_draws <- as_draws_array(fit_artL)
post_df    <- as_draws_df(fit_artL)

#summarise posterior draws
rhat_tbl <- posterior::summarise_draws(post_draws, mean, sd, ~ quantile(.x, c(0.025, 0.50, 0.975)),rhat, ess_bulk, ess_tail)
rhat_tbl

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

#WAIFW Beta Heatmaps
waifw_VT_fit_artL  <- summarise_matrix_param("beta_H_VT")
waifw_NVT_fit_artL <- summarise_matrix_param("beta_H_NVT")

#==================================================================== 
  
#import the posterior dataset
fit_artS <- rio::import(here::here("results", "fit_hhbayes_shortART.rds"))

#extract posterior draws
post_draws <- as_draws_array(fit_artS)
post_df    <- as_draws_df(fit_artS)

#summarise posterior draws
rhat_tbl <- posterior::summarise_draws(post_draws, mean, sd, ~ quantile(.x, c(0.025, 0.50, 0.975)),rhat, ess_bulk, ess_tail)
rhat_tbl

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

#WAIFW Beta Heatmaps
waifw_VT_fit_artS  <- summarise_matrix_param("beta_H_VT")
waifw_NVT_fit_artS <- summarise_matrix_param("beta_H_NVT")
  
#==================================================================== 

#import the posterior dataset
fit_neg <- rio::import(here::here("results", "fit_hhbayes_negHIV.rds"))

#extract posterior draws
post_draws <- as_draws_array(fit_neg)
post_df    <- as_draws_df(fit_neg)

#summarise posterior draws
rhat_tbl <- posterior::summarise_draws(post_draws, mean, sd, ~ quantile(.x, c(0.025, 0.50, 0.975)),rhat, ess_bulk, ess_tail)
rhat_tbl

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

#WAIFW Beta Heatmaps
waifw_VT_fit_neg  <- summarise_matrix_param("beta_H_VT")
waifw_NVT_fit_neg <- summarise_matrix_param("beta_H_NVT")


#==================================================================== 

sens_betaVT <- print((plot_waifw_5x5(waifw_VT_fit_main,  "A") + labs(subtitle = "Household VT Transmission (WAIFW \u03b2) if unknown adult HIV status is the same as that of HH index") | 
                      plot_waifw_5x5(waifw_VT_fit_artL, "B") + labs(subtitle = "Household VT Transmission (WAIFW \u03b2) if unknown adult HIV status is categorized as prolonged ART use"))/
                     (plot_waifw_5x5(waifw_VT_fit_artS,  "C") + labs(subtitle = "Household VT Transmission (WAIFW \u03b2) if unknown adult HIV status is categorized as shorter ART use") | 
                      plot_waifw_5x5(waifw_VT_fit_neg, "D") + labs(subtitle = "Household VT Transmission (WAIFW \u03b2) if unknown adult HIV status is categorized as negative"))) 

sens_betaNVT <- print((plot_waifw_5x5(waifw_NVT_fit_main,  "A") + labs(subtitle = "Household NVT Transmission (WAIFW \u03b2) if unknown adult HIV status is the same as that of HH index") | 
                      plot_waifw_5x5(waifw_NVT_fit_artL, "B") + labs(subtitle = "Household NVT Transmission (WAIFW \u03b2) if unknown adult HIV status is categorized as prolonged ART use"))/
                     (plot_waifw_5x5(waifw_NVT_fit_artS,  "C") + labs(subtitle = "Household NVT Transmission (WAIFW \u03b2) if unknown adult HIV status is categorized as shorter ART use") | 
                      plot_waifw_5x5(waifw_NVT_fit_neg, "D") + labs(subtitle = "Household NVT Transmission (WAIFW \u03b2) if unknown adult HIV status is categorized as negative"))) 

#combined plots
ggsave(here::here("output", "suppl_sensVT_WAIFW.png"), 
       plot = sens_betaVT, 
       width = 15, 
       height = 10, 
       unit = "in", 
       dpi = 300)

#combined plots
ggsave(here::here("output", "suppl_sensNVT_WAIFW.png"), 
       plot = sens_betaNVT, 
       width = 15, 
       height = 10, 
       unit = "in", 
       dpi = 300)
