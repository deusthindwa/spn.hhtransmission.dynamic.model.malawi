#written by Deus
#13/02/2024
#household pneumococcal carriage transmission in Malawi

#====================================================================

#load packages
pacman::p_load(char = c("lubridate", "tidyverse", "dplyr","msm", "here", "rio", "scales", "boot", "magrittr",  "mvtnorm", "zoo", "patchwork", "ggplotify", "sf",
                        "PropCIs", "reshape2","purrr", "msm", "minqa", "ggridges", "timetk", "ggbreak", "ggpubr", "gridExtra", "doParallel", "igraph", "rgdal"))

#====================================================================

#set seed for entire session globally to ensure reproducibility using a task call
addTaskCallback(function(...) {set.seed(12345);TRUE})

#turn off the global task call for set seed if needed
#removeTaskCallback(1)

#getting and putting datasets in right order for analysis
source(here("script", "1_data_wrangling.R"))
