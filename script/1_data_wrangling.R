#written by Deus+Ana
#23/01/2025
#household pneumococcal carriage transmission in Malawi

#====================================================================

#import spn overall household data
spn_hh <- 
  rio::import(here("data", "nasomuneHH_data.xlsx")) %>%
  dplyr::mutate(hiv = factor(if_else(hiv==0, 'hiv-',
                                        if_else(hiv==1, 'hiv+_artS',
                                                if_else(hiv==2, 'hiv+_artL', NA_character_)))),
                sex = factor(if_else(sex=='Female', 'female',
                                              if_else(sex=='Male', 'male', NA_character_))),
                agecat = factor(if_else(is.na(age), NA_character_,
                                        if_else(age <5, 'child', 'adult'))), #later make it 0-4,5-17, 18+
                pns = if_else(pns == 'neg', 1L, 2L),
                stg = if_else(stg == 'neg', 1L,
                              if_else(stg == 'VT', 2L, 3L)),
                
                denscat = if_else(is.na(st), 'none',
                                 if_else(!is.na(st) & dens < 26800, 'low', 'high'))) %>% #26800 median carriage density
  dplyr::arrange(pid, vno) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(vno = as.integer(if_else(!is.na(hiv), 1:n(), vno))) %>%
  dplyr::ungroup()

#====================================================================

#obtain spn baseline household data
spn_hhB <-
  dplyr::left_join(
    spn_hh %>%
      dplyr::select(hhid, pid, vno, pns, stg, age, agecat, sex, hiv, dens, denscat) %>%
      dplyr::filter(vno == 1),
    
    spn_hh %>%
      dplyr::select(hhid, pid, vno, pns, stg, age, agecat, sex, hiv, dens, denscat) %>%
      dplyr::filter(vno == 1) %>%
      dplyr::group_by(hhid) %>%
      tally() %>%
      dplyr::ungroup() %>%
      dplyr::rename('hhsize'='n'))

#====================================================================

#import spn patient level baseline and follow up data
spn_hhF <- 
  dplyr::left_join(
    spn_hh %>%
      dplyr::select(hhid, pid, vno, pns, stg, agecat, sex, hiv, dens, denscat),
    
    spn_hh %>%
      dplyr::select(hhid, pid, vno, pns, stg, agecat, sex, hiv, dens, denscat) %>%
      dplyr::filter(vno == 1) %>%
      dplyr::group_by(hhid) %>%
      tally() %>%
      dplyr::ungroup() %>%
      dplyr::rename('hhsize'='n'))
