#written by Deus
#10/04/2026
#WAIFW household pneumococcal carriage transmission modelling in Malawi

#====================================================================

#number of visits and carriage
spn_hhF %>%
  dplyr::group_by(vno, stg) %>%
  dplyr::tally() %>%
  dplyr::mutate(prev = n/sum(n), N = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(stg = if_else(stg==1, 'Uncol', if_else(stg==2, 'VT', 'NVT'))) %>%
  dplyr::ungroup()

#age distribution (mean/sd)
spn_hhB %>%
  dplyr::filter(!is.na(age)) %>%
  dplyr::mutate(stg = if_else(stg == 1, 'Uncol', if_else(stg == 2, 'VT', 'NVT'))) %>%
  dplyr::group_by(stg) %>%
  summarise(mX = mean(age), dX = sd(age))

#age grouped
spn_hhB %>%
  dplyr::group_by(agecatx, stg) %>%
  dplyr::tally() %>%
  dplyr::mutate(prev = n/sum(n),
                N = sum(n),
                stg = if_else(stg==1, 'Uncol', if_else(stg==2, 'VT', 'NVT'))) %>%
  dplyr::ungroup()

#hiv status
spn_hhB %>%
  dplyr::filter(agecatx == 'adult') %>%
  dplyr::group_by(hiv, stg) %>%
  tally() %>%
  dplyr::mutate(prev = n/sum(n), N = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(stg = if_else(stg==1, 'Uncol', if_else(stg==2, 'PCV13 serotypes', 'Non-PCV13 serotypes'))) %>%
  dplyr::ungroup()

#sex and carriage
spn_hhB %>%
  dplyr::filter(!is.na(stg)) %>%
  dplyr::group_by(sex, stg) %>%
  dplyr::tally() %>%
  dplyr::mutate(prev = n/sum(n), 
                N = sum(n),
                stg = if_else(stg==1, 'Uncol', if_else(stg==2, 'VT', 'NVT'))) %>%
  dplyr::ungroup()

#household size mean/sd
spn_hhB %>%
  dplyr::mutate(hhsizex = if_else(hhsize >=3 & hhsize <=5, '3-5',
                                  if_else(hhsize >=6 & hhsize <=20, '6+', '2'))) %>%
  dplyr::group_by(stg) %>%
  summarise(mX = mean(hhsize), dX = sd(hhsize)) %>%
  dplyr::ungroup()

# serotype group
spn_hhB %>%
  dplyr::mutate(hhsizex = if_else(hhsize >=3 & hhsize <=5, '3-5',
                                  if_else(hhsize >=6 & hhsize <=20, '6+', '2'))) %>%
  dplyr::group_by(hhsizex, stg) %>%
  dplyr::tally() %>%
  dplyr::mutate(prev = n/sum(n), N = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(stg = if_else(stg==1, 'Uncol', if_else(stg==2, 'PCV13 serotypes', 'Non-PCV13 serotypes'))) %>%
  dplyr::ungroup()

#household size
spn_hhB %>%
  dplyr::group_by(hhsize, stg) %>%
  dplyr::tally() %>%
  dplyr::mutate(prev = n/sum(n), N = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(stg = if_else(stg==1, 'Uncol', if_else(stg==2, 'PCV13 serotypes', 'Non-PCV13 serotypes'))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(stg != 'Uncol') 
  