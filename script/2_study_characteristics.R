#written by Deus+Ana
#23/01/2025
#household pneumococcal carriage transmission in Malawi

#====================================================================


#number of visits and carriage
vno <-
  spn_hhF %>%
  dplyr::group_by(vno, stg) %>%
  dplyr::tally() %>%
  dplyr::mutate(prev = n/sum(n), N = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(stg = if_else(stg==1, 'Uncol', if_else(stg==3, 'VT', 'NVT'))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(stg != 'Uncol') %>%
  
  ggplot(mapping = aes(x = vno, y = prev, color = stg, fill = stg)) + 
  geom_bar(stat = "identity", size = 0.7, position = position_stack(vjust = 0.5), color = "black") +
  geom_text(aes(label = n, fontface = 2), size = 5, color = "black", position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('#fc9272','#addd8e')) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(a)", x = "Visit number", y = "Carriage prevalence") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 0.60), labels = scales::percent_format(accuracy = 1)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))


#age distribution
agecat <- 
  spn_hhB %>%
  dplyr::filter(!is.na(age)) %>%
  dplyr::group_by(agecat, stg) %>%
  dplyr::tally() %>%
  dplyr::mutate(prev = n/sum(n),
                N = sum(n),
                stg = if_else(stg==1, 'Uncol', if_else(stg==3, 'VT', 'NVT'))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(stg != 'Uncol') %>%
  
  ggplot(mapping = aes(x = prev, y = agecat, color = stg, fill = stg)) + 
  geom_bar(stat = "identity", color = "black", size = 0.7) +
  geom_text(aes(label = n, fontface = 2), size = 4, color = "black", position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('#fc9272','#addd8e')) +
  theme_classic(base_size = 10, base_family = "American Typewriter") +
  labs(title = "", x = "Carriage prevalence", y = "Age group") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(legend.position = "none")

agedist <- 
  spn_hhB %>%
  dplyr::filter(!is.na(age)) %>%
  dplyr::mutate(stg = if_else(stg == 1, 'Uncol', if_else(stg == 3, 'VT', 'NVT'))) %>%
  dplyr::filter(stg != 'Uncol') %>%
  
  ggplot(aes(x = age, fill = stg, group = stg)) + 
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c('#d95f0e','#31a354')) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(b)", x = "Age distribution of participants", y = "Probability density") +
  scale_y_continuous(breaks = seq(0, 0.04, 0.01), limits = c(0, 0.04), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(0, 60, 10)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11)) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))


#sex and carriage
sex <-
  spn_hhB %>%
  dplyr::filter(!is.na(stg)) %>%
  dplyr::group_by(sex, stg) %>%
  dplyr::tally() %>%
  dplyr::mutate(prev = n/sum(n), 
                N = sum(n),
                stg = if_else(stg==1, 'Uncol', if_else(stg==3, 'VT', 'NVT'))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(stg != 'Uncol') %>%
  
  ggplot(mapping = aes(x = sex, y = prev, color = stg, fill = stg)) + 
  geom_bar(stat = "identity", color = "black", size = 0.7) +
  geom_text(aes(label = n, fontface = 2), size = 5, color = "black", position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('#fc9272','#addd8e')) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(c)", x = "Sex", y = "Carriage prevalence") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))


#density of spn carriage
denscatDS <-
  spn_hhF %>%
  dplyr::filter(!is.na(dens)) %>%
  dplyr::group_by(vno, denscat, stg) %>%
  dplyr::tally() %>%
  dplyr::mutate(prev = n/sum(n),
                N = sum(n),
                stg = if_else(stg==1, 'Uncol', if_else(stg==3, 'VT', 'NVT'))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(stg != 'Uncol')

densc1 <-
  denscatDS %>% dplyr::filter(vno == 1) %>%
  ggplot(mapping = aes(x = prev, y = denscat, color = stg, fill = stg)) + 
  geom_bar(stat = "identity", color = "black", size = 0.7) +
  geom_text(aes(label = n, fontface = 2), size = 4, color = "black", position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('#fc9272','#addd8e')) +
  theme_classic(base_size = 10, base_family = "American Typewriter") +
  labs(title = "", x = "Share of carriage", y = "Carriage density") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(legend.position = "none")

densc2 <-
  denscatDS %>% dplyr::filter(vno == 2) %>%
  ggplot(mapping = aes(x = prev, y = denscat, color = stg, fill = stg)) + 
  geom_bar(stat = "identity", color = "black", size = 0.7) +
  geom_text(aes(label = n, fontface = 2), size = 4, color = "black", position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('#fc9272','#addd8e')) +
  theme_classic(base_size = 10, base_family = "American Typewriter") +
  labs(title = "", x = "Share of carriage", y = "Carriage density") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(legend.position = "none")

densc3 <-
  denscatDS %>% dplyr::filter(vno == 3) %>%
  ggplot(mapping = aes(x = prev, y = denscat, color = stg, fill = stg)) + 
  geom_bar(stat = "identity", color = "black", size = 0.7) +
  geom_text(aes(label = n, fontface = 2), size = 4, color = "black", position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('#fc9272','#addd8e')) +
  theme_classic(base_size = 10, base_family = "American Typewriter") +
  labs(title = "", x = "Share of carriage", y = "Carriage density") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(legend.position = "none")

densc <- 
  spn_hhF %>%
  dplyr::filter(!is.na(dens)) %>%
  dplyr::mutate(stg = if_else(stg == 1, 'Uncol', if_else(stg == 3, 'VT', 'NVT')),
                vno = if_else(vno == 1, 'Visit 1', if_else(vno == 2, 'Visit 2', 'Visit 3'))) %>%
  
  ggplot(aes(x = dens, fill = stg, group = stg)) + 
  geom_density(alpha = 0.7) +
  scale_x_continuous(trans = log10_trans()) +
  scale_fill_manual(values = c('#d95f0e','#31a354')) +
  facet_grid(.~vno) + 
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(d)", x = "Carriage density (log_CFU/ml)", y = "Probability density") +
  scale_y_continuous(breaks = seq(0, 0.3, 0.05), limits = c(0,0.5), labels = scales::percent_format(accuracy = 1)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11)) + 
  theme(strip.background = element_rect(fill = "gray80"), strip.text.x = element_text(size=18)) +
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))


#household size and carriage
hhsize <- 
  spn_hhB %>%
  dplyr::group_by(hhsize, stg) %>%
  dplyr::tally() %>%
  dplyr::mutate(prev = n/sum(n), N = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(stg = if_else(stg==1, 'Uncol', if_else(stg==3, 'PCV13 serotypes', 'Non-PCV13 serotypes'))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(stg != 'Uncol') %>%
  
  
  ggplot(mapping = aes(x = hhsize, y = prev, color = stg, fill = stg)) + 
  geom_bar(stat = "identity", size = 0.7, position = position_stack(vjust = 0.5), color = "black") +
  geom_text(aes(label = n, fontface = 2), size = 5, color = "black", position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('#fc9272','#addd8e')) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(e)", x = "Number of household members", y = "Carriage prevalence") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(2, 9, 1)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  guides(fill = guide_legend(title = "Serotype group")) +
  theme(legend.position = 'right', legend.text = element_text(size = 16), legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"), legend.key.spacing.y = unit(0.2, "cm")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))


#combined plots
ggsave(here("output", "fig1_carr_char.png"),
       plot = ((vno | (agedist | inset_element(agecat, left = 0.3, right = 0.9,  bottom = 0.65, top = 0.95)) | sex | plot_layout(ncol = 3, width = c(2,3,1))) / 
              (densc | inset_element(densc1, left = 0.01, right = 0.3,  bottom = 0.65, top = 0.95) | inset_element(densc2, left = 0.35, right = 0.65, bottom = 0.65, top = 0.95) |
              inset_element(densc3, left = 0.70, right = 0.99, bottom = 0.65, top = 0.95)) | hhsize | plot_layout(ncol = 2, width = c(2,1))), 
       width = 24, height = 14, unit="in", dpi = 300)

