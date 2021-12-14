# Code to analyze meta data
# Em Lim
# Dec 2021, updated Oct 2026

# run renv::restore() if this is your first time spinning up this project
# renv::restore()

# Load packages
library(tidyverse)
library(ggplot2)
library(metafor)
library(renv)
library(emmeans) # for predictions
library(jtools) # for quick coeff plots
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(patchwork)

# Packages for map
library(sf)
library(rnaturalearth) # another simple world map?

# Source function file
source("Code/Functions.R")

# Set plotting theme
theme_set(theme_classic())

# Set preferences so I can see all my columns
rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

# List of papers generated in "Lit Search Script.R" then I added the meta-data to the following .csv

# Load data -------------------------------------------------------
# And do a million manipulations......
meta_full <- read_csv("Data/meta_analysis_data.csv") %>%
  rename_factors() %>% # function to rename pp, consumer, and nat bin factors
  mutate(
    # standard dev calcs
    sd_treatment = replace_na(sd_treatment, 999),# assigning a dummy to replace
    sd_control = replace_na(sd_control, 999),
    sd_treatment = ifelse(
      sd_treatment == 999,
      se_treatment * sqrt(treatment_n),
      sd_treatment),# replacing that dummy with the estimated SD
    sd_control = ifelse(sd_control == 999, se_control * sqrt(control_n), sd_control)# replacing that dummy with the estimated SE
    ) %>%
  filter(response != "protein") %>% # removes 3 rows
  filter(!(response == "count" & measured == "skeleton")) %>% # removes 1 row of skeletal density, i rename density = count with the rename_factors function
  # add sst data
  left_join(read_csv("Data/Output_data/meta_sst_final.csv") %>% rename(sst = annual_sst), by = "indiv") %>%
  # calculate added N
  rowwise() %>%
  mutate(
    # first calculate total N normally
    sum_all_N_treat = sum(nh4_treatment, no3_treatment, total_nitrogen_treatment, 
                          na.rm = TRUE),
    # then calculate total nh4 and no3
    sum_nh4_no3_treat = sum(nh4_treatment, no3_treatment, na.rm = TRUE),
    # then recalc total for rows where total = nh4 + no3
    sum_all_N_treat = 
      case_when(sum_nh4_no3_treat == total_nitrogen_treatment ~ total_nitrogen_treatment, 
                TRUE ~ as.numeric(sum_all_N_treat)),
    # now do the same for the control N
    sum_all_N_cont = sum(nh4_control, no3_control, total_nitrogen_control, 
                          na.rm = TRUE),
    # then calculate total nh4 and no3
    sum_nh4_no3_cont = sum(nh4_control, no3_control, na.rm = TRUE),
    # then recalc total for rows where total = nh4 + no3
    sum_all_N_cont = 
      case_when(sum_nh4_no3_cont == total_nitrogen_control ~ total_nitrogen_control, 
                TRUE ~ as.numeric(sum_all_N_cont)),
    added_N = na_if((sum_all_N_treat - sum_all_N_cont), 0),
    added_N = ifelse(prov_units == "umol/L" | 
                               prov_units == "umol/L/day" |
                               prov_units == "umol/L/h", added_N, NA),
    added_nh4 = na_if((sum(nh4_treatment, na.rm = TRUE) - 
                              sum(nh4_control, na.rm = TRUE)), 0),
    sum_all_N_cont_adj = ifelse(
      sum_all_N_cont == 0 & sum_all_N_treat != 0, 0.001, sum_all_N_cont),
    added_ratio = ifelse(sum_all_N_cont == 0 & sum_all_N_treat == 0, NA,
                         sum_all_N_treat/sum_all_N_cont_adj),
    added_percent = ifelse(sum_all_N_cont == 0 & sum_all_N_treat == 0, NA,
                           added_N/sum_all_N_cont_adj),
    added_N_log = ifelse(added_N < 0, NA, log(added_N)),
    po4 = as.numeric(po4),
    added_po4 = po4_treat - po4_control,
    added_po4 = ifelse(prov_units == "umol/L" | 
                       prov_units == "umol/L/day" |
                       prov_units == "umol/L/h", added_po4, NA),
    added_po4_log = ifelse(added_po4 <= 0, NA, log(added_po4)) 
         ) %>%
  ungroup() %>%
  # fix the Pfister paper, the cover slip experiment was microalgae, not phyto
  mutate(pp = case_when(paper_num == 33 ~ "Algae", TRUE ~ as.character(pp)),
         pp = as.factor(pp),
         year = ifelse(title == "Does Environmental Context Alter Bivalve and Snail Effects on Eelgrass Growth and Survival?", 2024, year)) %>%
  # remove mangroves + cordgrass, not enough vascular plants to look at this pp group
  filter(pp != "Plants") %>%
  # join the large marine area .csv produced by LME_calc.R
  left_join(read_csv("Data/Output_data/LME_df.csv"), by = "paper_num")


# Calculate effect sizes -------------------------------------------------------
meta_effects <- meta_full %>%
  escalc(measure = "ROM", data = meta_full, 
         m1i = mean_treatment, sd1i = sd_treatment, # log(treatment/control)
         m2i = mean_control, sd2i = sd_control,
         n1i = treatment_n, n2i = control_n) %>% 
  drop_na(yi) %>%
  filter(indiv != 245) # exclude huge SD outlier from Zhang paper
# throws an error because of log(-) or log(0) values, which we filter out
  # there's 9 cases where the mean control or treatment is 0
  # and 8 cases where the ratio is negative
  # so 17 NA's total
  # but only 9 cases where og effect sizes = NA because of a 0 (no og rows were negative)


# Dataframe of both growth and endpoint data, but without double rows for studies where I recorded both endpoint and growth data
# Basically, the "og" rows are the "correct" data type for the measured metric
# so if i could get growth, but i could also get an endpoint, then the "og" would be the growth row
# Unless the correct data form was just individual timepoints (eg. for chl a), in which case the time points would be "og" and the "extra" would be the after - before
meta_data <- meta_effects %>%
  filter(og_extra == "og")  %>% 
  scale_vars() %>% # call function created in Functions.R
  mutate(growth_rate = ifelse(growth_rate == "B", "Y", growth_rate)) # so I can check for growth rate effect

# Dataframe of just growth data
meta_grow <- meta_effects %>%
  filter(growth_rate != "N") %>%
  scale_vars()

# Dataframe of just end point data
meta_endpoint <- meta_effects %>%
  filter(growth_rate != "Y") %>%
  scale_vars()

# set df for mods
meta_df <- meta_data %>%
  droplevels() 
# write_csv(meta_data, "Data/Output_data/meta_df.csv")

# df with just rows with added N info
meta_df_a <- meta_data %>%
  drop_na(added_N) %>%
  droplevels()


# create reduced df to compare aic values
meta_df_aic <- meta_data %>%
  filter(bn_proxy == "measure") %>% 
  drop_na(duration_log)%>% 
  drop_na(added_N) %>%
  filter(added_N >= 0)%>%
  droplevels() 


# Multilevel models ------------------------------------------------------------
# Link to book with help:
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/metareg.html

# There is a special way through which we can tell rma.mv to assume such nested random effects. We do this using a slash (/) to separate the higher- and lower-level grouping variable. To the left of /, we put in the level 3 (cluster) variable. To the right, we insert the lower-order variable nested in the larger cluster. Therefore, the general structure of the formula looks like this: ~ 1 | cluster/effects_within_cluster

# Recommended to set test = “t” to use a test similar to the Knapp-Hartung method
# Recommended to set method = “REML” for restricted maximum likelihood
# Make sure categorical variables are factors

# Model output #
# Variance Components
  # random-effects variances calculated for each level of our model
  # sigma^2.1 shows the level 3 between-cluster variance (the between-study heterogeneity variance tao^2 in a conventional meta-analysis)
  # sigma^2.2 shows the variance within clusters (level 2)
  # nlvls is the number of levels within each level. X studies which contained X effect sizes

# Test for Residual Heterogeneity:
# If the test for residual heterogeneity is significant (p < 0.05), this possibly indicates that other moderators not considered in the model are influencing the vaccine effectiveness.

# Look at Test of Moderators to determine whether there are significant differences between the subgroups
# If this p-value < 0.05 then we can reject null hypothesis H0 : β1 = β2 = 0 based on the omnibus test (this means our moderators are better than a null model)

# The Model Results are printed within a meta-regression framework. This means that we cannot directly extract the estimates in order to obtain the pooled effect sizes within subgroups.


# Null model
mod_null <- rma.mv(yi = yi,
                   V = vi,
                   slab = author,
                   data = meta_df,
                   method = "REML",
                   mods = ~ 1,
                   random = ~ 1 | paper_num/indiv, 
                   test = "t")

summary(mod_null)

# Caluclate I^2
calc_I2(mod_null) # use function

# prep model for emmeans
prep_null <- emmeans_prep(mod_null) # use function to prep model for emmeans

# make table
pred_null_table <- as.data.frame(emmeans(prep_null, ~1, type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>%
  cbind((meta_df %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n)))) %>% 
  ungroup() %>%
  mutate(pp = "All")


### Latitude -------------------------------------------------------------------
# model with continuous latitude
mod_lat <- rma.mv(yi = yi,
                  V = vi,
                  slab = author,
                  data = meta_df,
                  method = "REML",
                  mods = ~ lat_stand*pp,
                  random = ~ 1 | paper_num/indiv, 
                  test = "t")
summary(mod_lat)
# write_csv(save_mod_output(mod_lat), "Tables/mod_lat.csv") 

# can we calculate the confidence intervals for the slopes using emtrends?
rg_l <- save_slopes_fun(
  (qdrg(object = mod_lat, 
        data = meta_df, 
        at = list(lat_stand = c(-.5, .5)
                  ))), "lat_stand"
  )
# write_csv(rg_l, "Tables/mod_lat_slopes.csv") 


# latitude categories
mod_tt <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = meta_df,
                 method = "REML",
                 mods = ~ tt*pp,
                 random = ~ 1 | paper_num/indiv, 
                 test = "t")
summary(mod_tt)
# write_csv(save_mod_output(mod_tt), "Tables/mod_tt.csv") 

# prep model for emmeans
prep_tt <- emmeans_prep(mod_tt) # use function to prep model for emmeans

# pairwise: did pp respond stronger temp vs tropics?
tt_pairs <- summary(pairs(emmeans(prep_tt,  ~ tt*pp), simple = "tt") ) %>%
  tidy_pairwise_output()
tt_pairs
# write_csv(tt_pairs, "Tables/tt_pairs.csv")

# pairwise: which pp responded stronger than others within each tt group
pairs(emmeans(prep_tt,  ~ tt*pp), simple = "pp") 

# prep table of estimated means
pred_tt_table1 <- as.data.frame(emmeans(prep_tt, specs="tt", by = "pp", type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>% # numbers 1 to 10, could have any column name here instead of df
  left_join((meta_df %>%
               group_by(pp, tt) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("pp", "tt")) %>% 
  replace_na(list(n_lab = "n = 0")) %>%
  ungroup()


### Background nitrogen -------------------------------------------------------
# Model with continuous background nitrogen
meta_df_n <- meta_df %>% filter(bn_proxy == "measure")

mod_bn <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = meta_df_n,
                 method = "REML",
                 mods = ~ bn_stand*pp,
                 random = ~ 1 | paper_num/indiv, 
                 test = "t")
summary(mod_bn)
# write_csv(save_mod_output(mod_bn), "Tables/mod_bn.csv")

# can we calculate the confidence intervals for the slopes using emtrends?
rg_b <- save_slopes_fun(
  (qdrg(object = mod_bn, data = meta_df_n, at = list(bn_stand = c(-.5, .5)))), 
  "bn_stand")
# write_csv(rg_b, "Tables/mod_bn_slopes.csv") 


# background nutrients as a categorical variable
# just use the papers that described nutrient levels in the paper!
mod_bn_cat <- rma.mv(yi = yi,
                       V = vi,
                       slab = author,
                       data = meta_df %>% drop_na(bn_cat),
                       method = "REML",
                       mods = ~ bn_cat*pp,
                       random = ~ 1 | paper_num/indiv, 
                       test = "t")

summary(mod_bn_cat)
# write_csv(save_mod_output(mod_bn_cat), "Tables/mod_bn_cat.csv") 

# prep model for emmeans
prep_bn_cat <- emmeans_prep(mod_bn_cat)

# pairwise comparisons
bn_pairs <- summary(pairs(emmeans(prep_bn_cat,  ~ bn_cat*pp), simple = "bn_cat") ) %>%
  tidy_pairwise_output()
bn_pairs
# write_csv(bn_pairs, "Tables/bn_pairs.csv")

pairs(emmeans(prep_bn_cat, ~ bn_cat*pp), simple = "pp") 

# make table of means
pred_bn_table1 <- as.data.frame(emmeans(prep_bn_cat, specs="bn_cat", by = "pp", type="response", tran = NA)) %>% # emmeans to pull estimates
  mutate(y_num = 1:length(df)) %>%
  left_join((meta_df %>%
               group_by(pp, bn_cat) %>%
               summarise(n = n()) %>%
               mutate(n_lab = sub("^","n = ",n))), by = c("pp", "bn_cat")) %>%
  replace_na(list(n_lab = "n = 0")) %>%
  ungroup() 

# Model with just bn_cat, no pp
mod_bn_cat_r <- rma.mv(yi = yi,
                       V = vi,
                       slab = author,
                       data = meta_df,
                       method = "REML",
                       mods = ~ bn_cat,
                       random = ~ 1 | paper_num/indiv, 
                       test = "t")
summary(mod_bn_cat_r) # no difference
#write_csv(save_mod_output(mod_bn_cat_r), "Tables/mod_bn_cat_r.csv") 


### Background phosphate -------------------------------------------------------
nrow(meta_df %>% drop_na(po4)) # 77 effect sizes (30%)
nrow(meta_df %>% drop_na(po4)%>% count(title)) # 13 studies (30%)

# how many rows for each pp?
meta_df %>%
  drop_na(po4) %>%
  count(pp)

# how many unique papers/po4 values for each pp
po4_count <- meta_df %>%
  drop_na(po4) %>%
  count(paper_num, title, po4, pp)

# make df for po4
meta_df_p <- meta_df %>% 
  drop_na(po4) %>%
  droplevels()


# Model with continuous background phosphate
mod_po4_pp <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = meta_df_p,
                 method = "REML",
                 mods = ~ po4_stand*pp,
                 random = ~ 1 | paper_num/indiv, 
                 test = "t")
summary(mod_po4_pp)
# write_csv(save_mod_output(mod_po4_pp), "Tables/mod_po4_pp.csv")

# just po4, no pp covariate
mod_po4 <- rma.mv(yi = yi,
                     V = vi,
                     slab = author,
                     data = meta_df_p,
                     method = "REML",
                     mods = ~ po4_stand,
                     random = ~ 1 | paper_num/indiv, 
                     test = "t")
summary(mod_po4)
# write_csv(save_mod_output(mod_po4), "Tables/mod_po4.csv")

### Temperature -------------------------------------------------------
# Model with continuous annual SST temperature
mod_sst <- rma.mv(yi = yi,
                  V = vi,
                  slab = author,
                  data = meta_df,
                  method = "REML",
                  mods = ~ sst_stand*pp, #+ added_N_stand,
                  random = ~ 1 | paper_num/indiv, 
                  test = "t")
summary(mod_sst)
# write_csv(save_mod_output(mod_sst), "Tables/mod_sst.csv")

# can we calculate the confidence intervals for the slopes using emtrends?
rg_s <- save_slopes_fun(
  (qdrg(object = mod_sst, data = meta_df, at = list(sst_stand = c(-.5, .5)))),
  "sst_stand")
# write_csv(rg_s, "Tables/mod_sst_slopes.csv") 

# Model with categorical annual SST temperature
mod_sst_cat <- rma.mv(yi = yi,
                  V = vi,
                  slab = author,
                  data = meta_df,
                  method = "REML",
                  mods = ~ sst_cat*pp,
                  random = ~ 1 | paper_num/indiv, 
                  test = "t")
summary(mod_sst_cat)
# write_csv(save_mod_output(mod_sst_cat), "Tables/mod_sst_cat.csv") 

# prep model for emmeans
prep_sst_cat <- emmeans_prep(mod_sst_cat)

# pairwise comparisons
sst_pairs <- summary(pairs(emmeans(prep_sst_cat,  ~ sst_cat*pp), simple = "sst_cat") ) %>%
  tidy_pairwise_output()
sst_pairs
# write_csv(sst_pairs, "Tables/sst_pairs.csv")

pairs(emmeans(prep_sst_cat,  ~ sst_cat*pp), simple = "pp") 


# table of means
pred_sst_table <- as.data.frame(emmeans(prep_sst_cat, specs="sst_cat", by = "pp", type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>% # numbers 1 to 10, could have any column name here instead of df
  left_join((meta_df %>%
               group_by(pp, sst_cat) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("pp", "sst_cat")) %>% 
  replace_na(list(n_lab = "n = 0")) %>%
  ungroup()


### rerun main mods with added N ----------------------
# model with continuous latitude
mod_lat_a <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = meta_df_a,
                    method = "REML",
                    mods = ~ lat_stand*pp + added_N_stand,
                    random = ~ 1 | paper_num/indiv, 
                    test = "t")
summary(mod_lat_a)
# write_csv(save_mod_output(mod_lat_a), "Tables/mod_lat_a.csv") 

# Model with continuous background_nutrients
mod_bn_a <- rma.mv(yi = yi,
                   V = vi,
                   slab = author,
                   data = meta_df_a %>% filter(bn_proxy == "measure"),
                   method = "REML",
                   mods = ~ bn_stand*pp + added_N_stand,
                   random = ~ 1 | paper_num/indiv, 
                   test = "t")
summary(mod_bn_a)
# write_csv(save_mod_output(mod_bn_a), "Tables/mod_bn_a.csv")

# Model with continuous annual SST temperature
mod_sst_a <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = meta_df_a,
                    method = "REML",
                    mods = ~ sst_stand*pp + added_N_stand,
                    random = ~ 1 | paper_num/indiv, 
                    test = "t")
summary(mod_sst_a)
# write_csv(save_mod_output(mod_sst_a), "Tables/mod_sst_a.csv")


## Covariate models ----------------------------------------------------------
### Experimental design ----------------------------------------------------
mod_nat <- rma.mv(yi = yi,
                        V = vi,
                        slab = author,
                        data = meta_df,
                        method = "REML",
                        mods = ~ nat_bin2,
                        random = ~ 1 | paper_num/indiv, 
                        test = "t")
summary(mod_nat)

# prep model for emmeans
prep_nat <- emmeans_prep(mod_nat) # use function to prep model for emmeans

# pairwise
summary(pairs(emmeans(prep_nat,  ~ nat_bin2)))

# make table
pred_nat_table <- as.data.frame(emmeans(prep_nat, specs="nat_bin2", type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>%
  left_join((meta_df %>%
               group_by(nat_bin2) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("nat_bin2")) %>% 
  ungroup() %>%
  mutate(pp = "All")


# experimental design vs pp
mod_nat_pp <- rma.mv(yi = yi,
                  V = vi,
                  slab = author,
                  data = meta_df,
                  method = "REML",
                  mods = ~ nat_bin2*pp,
                  random = ~ 1 | paper_num/indiv, 
                  test = "t")
summary(mod_nat_pp)
# write_csv(save_mod_output(mod_nat_pp), "Tables/mod_nat_pp.csv") 

# prep model for emmeans
prep_nat_pp <- emmeans_prep(mod_nat_pp) # use function to prep model for emmeans

# pairwise comparisons
nat_pairs <- summary(pairs(emmeans(prep_nat_pp,  ~ nat_bin2*pp), simple = "nat_bin2") ) %>%
  tidy_pairwise_output()
nat_pairs
# write_csv(nat_pairs, "Tables/nat_pairs.csv")

# make table
pred_nat_pp_table <- as.data.frame(emmeans(prep_nat_pp, specs="pp", by = c("nat_bin2"), type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>%
  left_join((meta_df %>%
               group_by(pp, nat_bin2) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("pp", "nat_bin2")) %>% 
  replace_na(list(n_lab = "n = 0")) %>%
  ungroup()

# plot
nat_plot <- effect_plot(pred_nat_pp_table, nat_bin2, 1) # call function created in Functions.R
nat_plot

### Inst vs growth vs endpoint -------------------------------------------------
mod_inst_grow <- rma.mv(yi = yi,
                        V = vi,
                        slab = author,
                        data = meta_df,
                        method = "REML",
                        mods = ~ grow_time_end,
                        random = ~ 1 | paper_num/indiv, 
                        test = "t")
summary(mod_inst_grow) # no diffs
# write_csv(save_mod_output(mod_inst_grow), "Tables/mod_inst_grow.csv") # S1

# prep model for emmeans
prep_grow <- emmeans_prep(mod_inst_grow) # use function to prep model for emmeans

# pairwise
grow_pairs_s <- summary(pairs(emmeans(emmeans_prep(mod_inst_grow),  ~ grow_time_end), simple = "grow_time_end") ) %>%
  mutate(pp = 1) %>%
  tidy_pairwise_output() %>% select(-pp)
grow_pairs_s
# write_csv(grow_pairs_s, "Tables/grow_pairs_s.csv")

# make table
pred_inst_grow_table <- as.data.frame(emmeans(prep_grow, specs="grow_time_end", type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>%
  left_join((meta_df %>%
               group_by(grow_time_end) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("grow_time_end")) %>% 
  ungroup()%>%
  mutate(pp = "All")

# rerun model with pp
mod_inst_grow_pp <- rma.mv(yi = yi,
                        V = vi,
                        slab = author,
                        data = meta_df,
                        method = "REML",
                        mods = ~ grow_time_end*pp,
                        random = ~ 1 | paper_num/indiv, 
                        test = "t")
summary(mod_inst_grow_pp)
# write_csv(save_mod_output(mod_inst_grow_pp), "Tables/mod_inst_grow_pp.csv") 

# prep model for emmeans
prep_grow_pp <- emmeans_prep(mod_inst_grow_pp) # use function to prep model for emmeans

# pairwise: between data types
grow_pairs <- summary(pairs(emmeans(prep_grow_pp,  ~ grow_time_end*pp), simple = "grow_time_end") ) %>%
  tidy_pairwise_output()
grow_pairs
# write_csv(grow_pairs, "Tables/grow_pairs.csv")
# For Algae and Coral, endpoint and growth rates differ
# algae = higher growth, coral = higher endpoint

# make table
pred_inst_grow_pp_table <- as.data.frame(emmeans(prep_grow_pp, specs="pp", by = c("grow_time_end"), type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>%
  left_join((meta_df %>%
               group_by(pp, grow_time_end) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("pp", "grow_time_end")) %>% 
  replace_na(list(n_lab = "n = 0")) %>%
  ungroup()

# plot
metric_plot <- effect_plot(pred_inst_grow_pp_table, grow_time_end, 1)
metric_plot

### Response -------------
mod_res <- rma.mv(yi = yi,
                   V = vi,
                   slab = author,
                   data = meta_df,
                   method = "REML",
                   mods = ~ res,
                   random = ~ 1 | paper_num/indiv, 
                   test = "t")
summary(mod_res)

# prep model for emmeans
prep_res <- emmeans_prep(mod_res) # use function to prep model for emmeans

# pairwise
summary(pairs(emmeans(prep_res,  ~ res )))

# make table
pred_res_table <- as.data.frame(emmeans(prep_res, specs="res", type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>%
  left_join((meta_df %>%
               group_by(res) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("res")) %>% 
  ungroup()%>%
  mutate(pp = "All")

# add pp to productivity response model
mod_res_pp <- rma.mv(yi = yi,
                      V = vi,
                      slab = author,
                      data = meta_df,
                      method = "REML",
                      mods = ~ res*pp,
                      random = ~ 1 | paper_num/indiv, 
                      test = "t")
summary(mod_res_pp)
# write_csv(save_mod_output(mod_res_pp), "Tables/mod_res_pp.csv") 

# prep model for emmeans
prep_res_pp <- emmeans_prep(mod_res_pp) # use function to prep model for emmeans

# pairwise: between responses
res_pairs <- summary(pairs(emmeans(prep_res_pp,  ~ res*pp), simple = "res") ) %>%
  tidy_pairwise_output() %>%
  drop_na(estimate)
res_pairs
# write_csv(res_pairs, "Tables/res_pairs.csv")

# make table
pred_res_pp_table <- as.data.frame(emmeans(prep_res_pp, specs="pp", by = c("res"), type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>%
  left_join((meta_df %>%
               group_by(pp, res) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("pp", "res")) %>% 
  replace_na(list(n_lab = "n = 0")) %>%
  ungroup()

# plot
res_plot <- effect_plot(pred_res_pp_table, res, 1)
res_plot

### Consumer taxa ----------------------------------------------------------
mod_taxa <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = meta_df,
                    method = "REML",
                    mods = ~ c_group,
                    random = ~ 1 | paper_num/indiv, 
                    test = "t")

summary(mod_taxa)

# prep model for emmeans
prep_taxa <- emmeans_prep(mod_taxa) # use function to prep model for emmeans

# pairwise
summary(pairs(emmeans(prep_taxa,  ~ c_group)))

# make table
pred_taxa_table <- as.data.frame(emmeans(prep_taxa, specs="c_group", type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>%
  left_join((meta_df %>%
               group_by(c_group) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("c_group")) %>% 
  ungroup()%>%
  mutate(pp = "All")

# add pp to consumer taxa level model
mod_taxa_pp <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = meta_df,
                    method = "REML",
                    mods = ~ c_group*pp,
                    random = ~ 1 | paper_num/indiv, 
                    test = "t")
summary(mod_taxa_pp)
# write_csv(save_mod_output(mod_taxa_pp), "Tables/mod_taxa_pp.csv") 

# prep model for emmeans
prep_taxa_pp <- emmeans_prep(mod_taxa_pp) # use function to prep model for emmeans

# pairwise: between consumer taxa levels
taxa_pairs <- summary(pairs(emmeans(prep_taxa_pp,  ~ c_group*pp), simple = "c_group") ) %>%
  tidy_pairwise_output() %>%
  drop_na(estimate)
taxa_pairs
# write_csv(taxa_pairs, "Tables/taxa_pairs.csv")

# make table
pred_taxa_pp_table <- as.data.frame(emmeans(prep_taxa_pp, specs="pp", by = c("c_group"), type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>%
  left_join((meta_df %>%
               group_by(pp, c_group) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("pp", "c_group")) %>% 
  replace_na(list(n_lab = "n = 0")) %>%
  ungroup()

# plot
c_taxa_plot <- effect_plot(pred_taxa_pp_table, c_group, 1)
c_taxa_plot



### Just primary producer groupings -------------------------------------------
mod_pp_no_int <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = meta_df,
                 method = "REML",
                 mods = ~ -1 + pp,
                 random = ~ 1 | paper_num/indiv, 
                 test = "t")
summary(mod_pp_no_int)
# write_csv(save_mod_output(mod_pp_no_int), "Tables/mod_pp.csv")

mod_pp <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = meta_df,
                 method = "REML",
                 mods = ~ pp,
                 random = ~ 1 | paper_num/indiv, 
                 test = "t")
summary(mod_pp)


# prep model for emmeans
prep_pp <- emmeans_prep(mod_pp) # use function to prep model for emmeans

# pairwise: between pp taxa levels
pp_pairs <- summary(pairs(emmeans(prep_pp,  ~ pp), simple = "pp") ) %>%
  drop_na(estimate)
pp_pairs
# write_csv(pp_pairs, "Tables/pppairs.csv")

# make table
pred_pp_table <- as.data.frame(emmeans(prep_pp, specs="pp", type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df),
         dum = "Overall") %>%
  left_join((meta_df %>%
               group_by(pp) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("pp")) %>% 
  replace_na(list(n_lab = "n = 0")) %>%
  ungroup()

# plot
pp_plot <- effect_plot(pred_pp_table, dum, 1)
pp_plot


### Duration -------------------------------------------------------------------
mod_dur <- rma.mv(yi = yi,
                  V = vi,
                  slab = author,
                  data = meta_df,
                  method = "REML",
                  mods = ~ duration_log*pp,
                  random = ~ 1 | paper_num/indiv, 
                  test = "t")
summary(mod_dur)
# write_csv(save_mod_output(mod_dur), "Tables/mod_dur2.csv")

ggplot(data= meta_df, aes(duration_log, yi, colour = pp, group = pp)) +
  geom_point() +
  geom_smooth(method = lm)

# prep for emmeans
prep_dur <- emmeans_prep(mod_dur) # use function to prep model for emmeans


# DURATION GLMS
# Does duration vary with latitude???
ggplot(data= meta_df, aes(lat_center, duration_days, colour = pp)) +
  geom_point() +
  geom_smooth(method = lm)

# model duration vs lat
mod_lm_dur_lat <- glmmTMB(duration_days ~ lat_center + (1|paper_num), 
                          family = tweedie(link = "log"),
                          data = meta_df)
plot(simulateResiduals(mod_lm_dur_lat))
summary(mod_lm_dur_lat) # looks like duration decreases with latitude

# model duration vs pp
mod_lm_dur_pp <- glmmTMB(duration_days ~ pp + (1|paper_num), 
                          family = tweedie(link = "log"),
                          data = meta_df)
plot(simulateResiduals(mod_lm_dur_pp))
summary(mod_lm_dur_pp)

pairs(emmeans(mod_lm_dur_pp, ~ pp)) # but really coral studies are just way longer

# model duration vs lat vs pp
mod_lm_dur_lat_pp <- glmmTMB(duration_days ~ lat_center*pp + (1|paper_num), 
                         family = tweedie(link = "log"),
                         data = meta_df)
plot(simulateResiduals(mod_lm_dur_lat_pp))
summary(mod_lm_dur_lat_pp)
# ok yes confirmed, there's no relationship between duration and lat, it's actually just that low lat coral studies are super long

# save output
# save_glm_output(mod_lm_dur_lat_pp, "Tables/dur_lat_glm2.csv")


AIC(mod_lm_dur_lat, mod_lm_dur_pp, mod_lm_dur_lat_pp)
# when I compare these with AIC, the just pp mod is best

# does duration vary with nat bin?
ggplot(data= meta_df, aes(nat_bin2, duration_days)) +
  geom_point() +
  stat_summary(fun = "mean", size = 1, alpha = 0.7) +
  stat_summary(geom = "errorbar", size = 0.5, alpha = 0.5) +
  facet_wrap(~ pp, scales = "free_y")

# model
mod_lm_dur_nat <- glmmTMB(duration_days ~ nat_bin2 + (1|paper_num), 
                          family = tweedie(link = "log"),
                          data = meta_df )
plot(simulateResiduals(mod_lm_dur_nat))
summary(mod_lm_dur_nat) # not really a big diff

# model
mod_lm_dur_pp <- glmmTMB(duration_days ~ pp + (1|paper_num), 
                          family = tweedie(link = "log"),
                          data = meta_df )
plot(simulateResiduals(mod_lm_dur_pp))
summary(mod_lm_dur_pp) # coral longer than others

# model
mod_lm_dur_nat_pp <- glmmTMB(duration_days ~ nat_bin2*pp + (1|paper_num), 
                         family = tweedie(link = "log"),
                         data = meta_df )
plot(simulateResiduals(mod_lm_dur_nat_pp))
summary(mod_lm_dur_nat_pp)
# save output
# save_glm_output(mod_lm_dur_nat_pp, "Tables/dur_nat_glm2.csv")


# pairwise: between experimental design levels
nat_dur_pairs <- summary(pairs(emmeans(mod_lm_dur_nat_pp,  ~ nat_bin2*pp), simple = "nat_bin2") ) %>%
  tidy_glm_pairwise_output() %>%
  drop_na(estimate)
nat_dur_pairs
# write_csv(taxa_pairs, "Tables/taxa_pairs.csv")

AIC(mod_lm_dur_nat, mod_lm_dur_pp, mod_lm_dur_nat_pp) # nat and pp best

emmeans(mod_lm_dur_nat_pp, ~ pp*nat_bin2)
pairs(emmeans(mod_lm_dur_nat_pp, ~ nat_bin2*pp), simple = "nat_bin2") 


# duration vs consumer
mod_lm_dur_c <- glmmTMB(duration_days ~ c_group + (1|paper_num), 
                        family = tweedie(link = "log"),
                        data = meta_df )
summary(mod_lm_dur_c)
summary(pairs(emmeans(mod_lm_dur_c,  ~ c_group)) ) # no diffs


# tt? ignore this, the outliers leverage the mean so much
ggplot(data= meta_df, aes(tt, duration_days)) +
  geom_point() +
  stat_summary(fun = "mean", size = 1, alpha = 0.7) +
  stat_summary(geom = "errorbar", size = 0.5, alpha = 0.5) +
  facet_wrap(~ pp)

# what is the nature of these dur = 0 studies
dur <- meta_df %>%
  filter(duration_days == 0) %>%
  count(nat_bin2, inst_time) 


### Added nitrogen -------------------------------------------------------------
meta_added <- meta_df %>% 
  drop_na(added_N) %>%
  filter(added_N >= 0) 

meta_add_ratio <- meta_df %>% 
  drop_na(added_ratio) %>%
  filter(added_ratio >= 0)

# how many?
nrow(meta_added %>% count(paper_num)) # 20

# run model
mod_add <- rma.mv(
  yi = yi,
  V = vi,
  slab = author,
  data = meta_added,
  method = "REML",
  mods = ~ added_N_log*pp,
  random = ~ 1 | paper_num / indiv,
  test = "t")
summary(mod_add)
# write_csv(save_mod_output(mod_add), "Tables/mod_add.csv")

# can we calculate the confidence intervals for the slopes using emtrends?
rg_a <- save_slopes_fun(
  (qdrg(object = mod_add, data = meta_added, at = list(added_N_log = c(-.5, .5)))),
  "added_N_log")
# write_csv(rg_a, "Tables/mod_added_slopes.csv") 

# run model
mod_add2 <- rma.mv(
  yi = yi,
  V = vi,
  slab = author,
  data = meta_add_ratio %>% filter(bn_proxy == "measure"),
  method = "REML",
  mods = ~ added_ratio*pp + bn_stand*pp,
  random = ~ 1 | paper_num / indiv,
  test = "t")
summary(mod_add2)


# model added_N vs lat
mod_lm_add_lat <- glmmTMB(added_N ~ lat_center*pp + (1|paper_num), 
                          family = tweedie(link = "log"),
                          data = meta_added )
# plot(simulateResiduals(mod_lm_add_lat))
summary(mod_lm_add_lat)
# save_glm_output(mod_lm_add_lat, "Tables/add_lat_glm.csv")

pairs(emmeans(mod_lm_add_lat, ~ pp)) 

ggplot(meta_added, aes(lat_center, added_N, colour = pp)) +
  geom_point()+
  geom_smooth(method = lm)

# added N vs experimental design?
mod_lm_add_nat <- glmmTMB(added_N ~ nat_bin2*pp + (1|paper_num), 
                          family = tweedie(link = "log"),
                        data = meta_added)
#plot(DHARMa::simulateResiduals(mod_lm_add_nat))
summary(mod_lm_add_nat)
# save_glm_output(mod_lm_add_nat, "Tables/add_nat_glm.csv")


pairs(emmeans(mod_lm_add_nat, ~ nat_bin2*pp), simple = "nat_bin2") 

# Does added N vary between pp*tt
mod_lm_add_tt <- glmmTMB(added_N ~ pp*tt + (1|indiv/paper_num), 
                         #     family = tweedie(link = "log"),
                         data = meta_added )
# plot(simulateResiduals(mod_lm_add))
summary(mod_lm_add_tt)

pairs(emmeans(mod_lm_add_tt, ~ tt*pp), simple = "tt") 


# added N vs consumer
mod_lm_add_c <- glmmTMB(added_N ~ c_group + (1|paper_num), 
                          family = tweedie(link = "log"),
                          data = meta_added )
#plot(simulateResiduals(mod_lm_add_c))
summary(mod_lm_add_c)
# save_glm_output(mod_lm_add_c, "Tables/mod_lm_add_c.csv")

summary(pairs(emmeans(mod_lm_add_c,  ~ c_group)) ) # jellies pee less than everything



### Added phosphate -------------------------------------------------------------
meta_added_p <- meta_df %>% 
  drop_na(added_po4) %>%
  # phyto is the only pp with more than one study
 # filter(pp == "Phytoplankton") %>% 
  droplevels()

# how many?
nrow(meta_added_p %>% count(paper_num)) # 11 studies, 75 effect sizes
d <- meta_added_p %>% count(pp, added_po4)

# run model
mod_add_p <- rma.mv(
  yi = yi,
  V = vi,
  slab = author,
  data = meta_added_p,
  method = "REML",
  mods = ~ added_po4,
  random = ~ 1 | paper_num / indiv,
  test = "t")
summary(mod_add_p)


## Epiphyte vs seagrass -------------------------------------------------------

# make a df of just seagrass effect sizes
s_nums <- meta_df %>%
  filter(pp_taxa == "S")

# just the fleshy algae effect sizes
e_nums <- meta_df %>%
  filter(pp_taxa == "A")

# of the seagrass papers, which of those papers also recorded a fleshy algae 
epi_sea_nums <- s_nums %>%
  filter(paper_num %in% e_nums$paper_num)

# ok now make a final df of those papers
meta_sge <- meta_df %>%
  filter(paper_num %in% epi_sea_nums$paper_num) %>%
  droplevels() %>%
  as.data.frame() 
# I checked all of the papers on algae here are actually epiphytes
# Also checked no epiphytes were tagged as phytoplankton

# just compare groups
mod_sge <-  rma.mv(yi = yi,
                       V = vi,
                       slab = author,
                       data = meta_sge,
                       method = "REML",
                       mods = ~ -1 + pp,
                       random = ~ 1 | paper_num/indiv, 
                       test = "t")
summary(mod_sge)

# model climate group
mod_lat_sge <-  rma.mv(yi = yi,
                      V = vi,
                      slab = author,
                      data = meta_sge,
                      method = "REML",
                      mods = ~ lat_center*pp,
                      random = ~ 1 | paper_num/indiv, 
                      test = "t")
summary(mod_lat_sge)
# write_csv(save_mod_output(mod_lat_sge), "Tables/mod_lat_sge.csv") 

# can we calculate the confidence intervals for the slopes using emtrends?
rg_l_sge <- save_slopes_fun(
  (qdrg(object = mod_lat_sge, data = meta_sge, at = list(lat_center = c(-.5, .5)))),
  "lat_center")
# write_csv(rg_l_sge, "Tables/mod_lat_sge_slopes.csv") 


# now run the model for background nutrients
mod_bn_sge <-  rma.mv(yi = yi,
                      V = vi,
                      slab = author,
                      data = meta_sge,
                      method = "REML",
                      mods = ~ pp*bn_cat,
                      random = ~ 1 | paper_num/indiv, 
                      test = "t")

summary(mod_bn_sge)
# write_csv(save_mod_output(mod_bn_sge), "Tables/mod_bn_sge.csv") 

# prep model for emmeans
prep_sge_bn <- emmeans_prep(mod_bn_sge) 

# pairwise: did pp respond stronger low vs high nutrients?
bn_pairs_sge <- summary(pairs(emmeans(prep_sge_bn,  ~ bn_cat*pp), simple = "bn_cat") ) %>%
  tidy_pairwise_output()
bn_pairs_sge
# write_csv(bn_pairs_sge, "Tables/bn_pairs_sge.csv")

# pp rank within region
pairs(emmeans(prep_sge_bn,  ~ bn_cat*pp), simple = "pp") 


# now run the model for SST
mod_sst_sge <-  rma.mv(yi = yi,
                      V = vi,
                      slab = author,
                      data = meta_sge,
                      method = "REML",
                      mods = ~ sst_center*pp,
                      random = ~ 1 | paper_num/indiv, 
                      test = "t")

summary(mod_sst_sge)
# write_csv(save_mod_output(mod_sst_sge), "Tables/mod_sst_sge.csv") 

# can we calculate the confidence intervals for the slopes using emtrends?
rg_s_sge <- save_slopes_fun((
  qdrg(object = mod_sst_sge, data = meta_sge, at = list(sst_center = c(-.5, .5)))),
  "sst_center")
# write_csv(rg_s_sge, "Tables/mod_sst_sge_slopes.csv") 

## Region ----------------------------------------------------------------------
mod_lme <- rma.mv(yi = yi,
                  V = vi,
                  slab = author,
                  data = meta_df,
                  method = "REML",
                  mods = ~ lme_name,
                  random = ~ 1 | paper_num/indiv, 
                  test = "t")
summary(mod_lme)

# prep model for emmeans
prep_lme <- emmeans_prep(mod_lme) # use function to prep model for emmeans
lme_pairwise <- summary(pairs(emmeans(prep_lme,  ~ lme_name))) %>%
  dplyr::transmute(contrast = contrast,
    estimate = ifelse(abs(estimate) < 0.005, sprintf("%.3f",estimate), sprintf("%.2f",estimate)),
         SE = ifelse(abs(SE) < 0.005, sprintf("%.3f",SE), sprintf("%.2f",SE)),
         `t-ratio` = ifelse(abs(t.ratio) < 0.005, sprintf("%.3f",t.ratio), sprintf("%.2f",t.ratio)),
         `p-value` = ifelse(p.value < 0.001, sprintf("< 0.001"), 
                          ifelse(p.value < 0.005, sprintf("%.3f",p.value), sprintf("%.2f",p.value))))
# write_csv(lme_pairwise, "Tables/lme_pairs.csv")



# Figures ----------------------------------------------------------------------

# Possible figure sizes: single column = 79mm, 2/3rd column = 110mm, double column = 168mm, maximum height of figure = 230mm.
# Font size in figures should be 8 pt.
# save in .eps or .pdf format
# 600 dpi 

# Define palettes
blue <- paste("#b9d1df", sep="")

# create some palettes
pal <- viridis::viridis(6)
pie(rep(1, 6), col = pal)
pal_tt <- c(pal[2], pal[4])
pal5 <- pal[2:6]
pal_i <- viridis::magma(10)
pal1 <- pal_i[6]

# just use the same pal as the rest of the figs
pal_map <- viridis::viridis(4)

# Labels
nat_labs <- c("Unstocked", "Stocked", "Mesocosm")
pp_labs <- c("Algae", "Coral", "Phytoplankton", "Plants")

# Maybe some help here?
# https://www.metafor-project.org/doku.php/tips:meta_regression_with_log_rr

## Figure 1 -------------------------------------------------------------------
### Fig.1a Lat vs temp vs nutrients ---------------------------------------------
set.seed(123)

ggplot() +
  # add NA's for bn first
  geom_jitter(data = meta_df %>% 
                filter(bn_proxy == "proxy"), 
              aes(abs_lat, sst),
              size = 0.5, alpha = 0.5, height = 1, width = 2.5,
              colour = "darkgrey") +
  # then add measured values
  geom_jitter(data = meta_data %>% 
                filter(bn_proxy == "measure"), 
              aes(abs_lat, sst, colour = background_nut_umolL),
              size = 0.5, alpha = 0.7, height = 1, width = 2.5)  +
  viridis::scale_colour_viridis() +
  viridis::scale_fill_viridis() +
  labs(x = "Absolute latitude", colour = expression("N"~(mu*M)), y = expression("Annual SST " ( degree*C))) +
  xlim(c(0, 80)) +
  pub_theme()
# ggsave("Figures/Fig.1a.png", device = "png", height = 60, width = 84, units = c("mm"), dpi = 600)


### Fig.1b - study site map --------
sf_use_s2(FALSE)

# Good sites for plotting help
# https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
# https://mhallwor.github.io/_pages/Tidyverse_intro

lat_long <- meta_data %>%
  distinct(paper_num, lat_dec, long_dec, tt) %>%
  as.data.frame() %>%
  mutate(lat = lat_dec,
         long = long_dec) %>%
  st_as_sf(coords = c("long_dec", "lat_dec")) %>%
  st_set_crs(4326) %>%
  st_transform(crs = "+proj=moll")

# Transform to a crs that uses meters as the distance
coords_jitter <- lat_long %>%
  st_jitter(nc_points, amount = 150000)
# st_transform(4326)

# Load simple world map from package "rnaturalearth"
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = "+proj=moll")

# set seed
set.seed(987654321)

# Create the elliptical mask with updated parameters
globe_mask_ellipse <- make_ellipse()

# Create an sf object
globe_mask_sf <- st_sf(geometry = globe_mask_ellipse)

# Clip world and points to the new elliptical mask
world_clipped <- st_intersection(world, globe_mask_sf)
points_clipped <- st_intersection(coords_jitter, globe_mask_sf)

# Plot with Mollweide and correct ellipse proportions
ggplot() +
  geom_sf(data = globe_mask_sf, fill = blue, color = NA) +
  geom_sf(data = world_clipped, fill = "white", color = "black") +
  # LME polygons
 # geom_sf(data = lme_clipped, aes(fill = LME_NAME), color = "black", alpha = 0.5) +
  # add coordinates for each study
  geom_sf(data = points_clipped, pch = 21, colour = "black", fill = pal1, size = 2.5, alpha = 0.6) +
  # Add the black circle (ellipse) around the Earth
  geom_sf(data = globe_mask_ellipse, fill = NA, color = "black", size = 0.5) +  # Black circle outline
  # Circular labels for LMEs
#  geom_sf_label(
#    data = lme_c,
#    aes(label = LME_NUMBER),   # replace with your label column
#    fill = "white",
#    color = "black",
#    label.r = unit(0.5, "lines"),  # rounded corners
#    label.size = 0.1,
#    size = 1.1) +
  # Set the coordinate system to Mollweide and let it manage aspect ratio
  coord_sf(crs = st_crs("+proj=moll"), 
           expand = FALSE) +  # No manual scaling, just preserve aspect ratio
  pub_theme() +
  theme(
    axis.text = element_blank(),    # Remove axis labels
    axis.title = element_blank(),   # Remove axis titles
    axis.ticks = element_blank(),    # Remove axis ticks
    axis.line = element_blank(),    # Remove the axis lines
    panel.border = element_blank(),  # Remove the black border around the plot
    legend.position = "none"
  ) +
  scale_fill_manual(values = viridis::viridis(66))
# ggsave("Figures/Fig.1b.png", width = 7.08661, height = 3, dpi = 600)


### Fig.1c - trend in publishing over time ------
meta_full2 <- meta_full %>% # include all papers that passed criteria
  mutate(year = case_when(paper_num == "45" ~ 2024,
                          paper_num == "44" ~ 2019,
                          TRUE ~ as.numeric(year)))%>%
  select(paper_num, year, tt) %>%
  unique()

# make df of all years
years <- min(meta_full2$year):max(meta_full2$year)
n_years <- length(years)

# add up cumulative sums
meta_y <- data.frame(
  year = rep(years, 2)) %>%
  mutate(tt = rep(c("Tropical", "Temperate"), each = n_years)) %>%
  left_join(meta_full2 %>%
              count(tt, year),
            by = c("year", "tt")) %>%
  rbind(
    data.frame(
      year = rep(years, 1)) %>%
      mutate(tt = rep(c("Total"), each = n_years)) %>%
      left_join(meta_full2 %>%
                  count(year),
                by = "year")
  ) %>%
  mutate(n = tidyr::replace_na(n, 0)) %>%
  group_by(tt) %>%
  mutate(cum_sum = cumsum(n)) %>%
  ungroup() %>%
  mutate(tt = factor(tt, levels = c("Total", "Temperate", "Tropical"))) %>%
  filter(!(tt == "Temperate" & cum_sum == 0))

# Plot
ggplot(meta_y, aes(x = year, y = cum_sum, colour = tt, lty = tt)) +
  geom_line(linewidth = 0.4) +
  scale_colour_manual(
    values = c(
      "Total"      = "black",
      "Tropical"   = "#22A884FF",
      "Temperate"  = "#414487FF")) +
  labs(x = "Publication Year", y = "Number of Publications", colour = "", lty = "")  +
  pub_theme() +
  theme(legend.position = c(0.3, 0.75),
        legend.title = element_blank())

# ggsave("Figures/Fig.1c.png", device = "png", height = 60, width = 84, units = c("mm"), dpi = 600)


## Fig.2 - pp x nat x tt --------------------------------------------------------
meta_counts <- meta_df %>%
  count(c_group, pp, nat_bin2, tt) %>% 
  mutate(pp = case_when(pp == "Endosymbiont" ~ "Endo",
                        pp == "Phytoplankton" ~ "Phyto",
                        TRUE ~ as.factor(pp)),
         pp = as.factor(pp),
         pp_num = as.numeric(pp),
         # nudge specific points slightly upward (by 0.1)
         pp_adj = case_when(
           c_group == "Molluscs" & pp == "Seagrass" & tt == "Tropical"  ~ pp_num + 0.2,
           c_group == "Molluscs" & pp == "Seagrass" &
             tt == "Temperate" & nat_bin2 == "Stocked" ~ pp_num - 0.1,
           c_group == "Sea cucumbers" & pp == "Phyto" & tt == "Tropical" ~ pp_num - 0.1,
           c_group == "Molluscs" & pp == "Phyto" & 
             tt == "Temperate" & nat_bin2 == "Mesocosm" ~ pp_num - 0.1,
           c_group == "Molluscs" & pp == "Algae" & tt == "Tropical"  ~ pp_num + 0.05,
           TRUE ~ pp_num),
         c_num = as.numeric(c_group),
         c_group_adj = case_when(
           c_group == "Fishes" & pp == "Algae" ~ c_num + 0.2,
           c_group == "Sponges" & pp == "Algae" ~ c_num + 0.2,
           c_group == "Birds" & pp == "Algae" ~ c_num + 0.2,
           c_group == "Birds" & pp == "Coral" ~ c_num + 0.25,
           c_group == "Fishes" & pp == "Endo" & nat_bin2 == "Unstocked" ~ c_num - 0.2,
           c_group == "Fishes" & pp == "Coral" & nat_bin2 == "Unstocked" ~ c_num - 0.1,
           TRUE ~ c_num
         ))

# plot
seed <- 11
ggplot() +
  geom_point(data = meta_counts, 
             aes(c_group_adj, pp_adj, fill = tt, pch = nat_bin2, size = n),
              alpha = 0.4, 
              position = position_jitter(seed = seed, height = 0.3, width = 0.3)) +
  geom_text(data = meta_counts,
            aes(x = c_group_adj, y = pp_adj, label = n, group = tt), # change to 1 if i add tt
            color = "black", size = 3, 
            position = position_jitter(seed = seed, height = 0.3, width = 0.3)) +
  labs(x = "Consumer group", y = "Primary producer group",
       pch = "Experiment", fill = "Climate", size = "n") +
  scale_fill_manual(values = (pal_tt)) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  guides(fill = guide_legend(override.aes = 
                               list(size = 3.5, pch = 21, alpha = 0.4)),
         pch = guide_legend(override.aes = 
                              list(size = 3.5, colour = "black", alpha = 0.9))) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.spacing = unit(0.1, "cm"),
          legend.box.spacing = unit(0.1, "cm")) + 
    scale_size_continuous(range = c(3,8),
                          breaks = c(5, 15, 25)) +
  scale_y_continuous(
    breaks = 1:length(levels(meta_counts$pp)),
    labels = levels(meta_counts$pp))+
  scale_x_continuous(
    breaks = 1:length(levels(meta_counts$c_group)),
    labels = levels(meta_counts$c_group))

# ggsave("Figures/Fig.2_pp_x_c.png", device = "png", height = 100, width = 140, units = c("mm"), dpi = 600)




## Fig.3 Overall estimates for pp groups -------------------------------------
# make mega df of all of the mod effect sizes! will use this in fig.5 too
# then just filter for the "all" effect sizes for this plot

mega_df <- 
  # add overall null model
  pred_null_table %>%
  mutate(type = "All") %>%
  rename(variable = 1) %>%
  # add overall pp
  rbind(pred_pp_table %>%
          mutate(type = "All") %>%
          rename(variable = dum)) %>%
  # experimental design all
  rbind(pred_nat_table %>%
          mutate(type = "Design") %>%
          rename(variable = nat_bin2)) %>%
  # experimental design pp
  rbind(pred_nat_pp_table %>%
          mutate(type = "Design") %>%
          rename(variable = nat_bin2)) %>%
  # add inst_grow overall 
  rbind(pred_inst_grow_table %>%
          mutate(type = "Measured") %>%
          rename(variable = grow_time_end)) %>%
  # add inst_grow pp
  rbind(pred_inst_grow_pp_table %>%
          mutate(type = "Measured") %>%
          rename(variable = grow_time_end)) %>%
  # add response
  rbind(pred_res_table %>%
          mutate(type = "Response") %>%
          rename(variable = res)) %>%
  # add response pp
  rbind(pred_res_pp_table %>%
          mutate(type = "Response") %>%
          rename(variable = res)) %>%
  # add c_taxa overall
  rbind(pred_taxa_table %>%
          mutate(type = "Consumer") %>%
          rename(variable = c_group)) %>%
  # add c_taxa pp
  rbind(pred_taxa_pp_table %>%
          mutate(type = "Consumer") %>%
          rename(variable = c_group)) %>%
  # tidy variables
  mutate(variable = str_to_sentence(variable),
         variable = case_when(variable == "Timepoint" ~ "Instantaneous",
                              variable == "Growth_rate_data" ~ "Growth rate",
                              variable == "Endpoint" ~ "End point",
                              variable == "Chl_a" ~ "Chl a",
                              TRUE ~ as.character(variable)
         ),
         type = factor(as.factor(type), 
                       levels = c("All", "Design", "Measured", "Response", "Consumer")),
         
         # nudge text position
         text_pos = case_when(
           pp == "Phytoplankton" & 
             type == "Design" & 
             variable == "Unstocked" ~ emmean - 0.2,
           pp == "Seagrass" & 
             type == "Design" & 
             variable == "Mesocosm" ~ emmean - 0.5,
           pp == "Seagrass" & 
             type == "Response" & 
             variable == "Count" ~ emmean + 0.1,
           pp == "Phytoplankton" & 
             type == "Consumer" & 
             variable == "Birds" ~ emmean - 0.3,
           TRUE ~ emmean),
         # reorder pp levels
         pp = factor(pp, levels =
                       c("All", "Algae", "Coral", "Endosymbiont", "Phytoplankton", "Seagrass")),
         point_size = ifelse(type == "All" & pp == "All", "small", "normal")
  ) 


# just make the overall figure, no covariates
overall_mega_df <- mega_df %>% 
  filter(type == "All") %>%
  mutate(variable = ifelse(pp == "All", "All", "Primary producer taxa"),
         pp = ifelse(pp == "All", "Overall", as.character(pp)),
         pp = factor(pp, levels = c("Seagrass", "Phytoplankton", "Endosymbiont", "Coral", "Algae", "Overall")))

# plot
ggplot(overall_mega_df, aes(emmean, pp, size = point_size)) +
  geom_point() +
  scale_size_manual(values = c("small" = 1, "normal" = 1.3)) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  facet_grid(
    rows = vars(variable),
    scales = "free_y",
    space = "free_y") +
  pub_theme() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)") +
  #scale_color_manual(values = (pal)) +
  # Add text labels inside the plots
  geom_text(aes(x = text_pos, label = n), # change to 1 if i add tt
            color = "black", size = 3, vjust = -0.5, hjust = 0)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

# ggsave("Figures/Fig.3.pp_overall.png", device = "png", height = 80, width = 90, units = c("mm"), dpi = 600)



## Fig.4 - Cont mods ----------------------------------------------------------

### Continuous latitude --------------------------------------------------------
slope_lat <- coef(summary(mod_lat))[2,1] # extract the slope

# Calculate slope for each pp
int_effect_sum_l <- slope_coeff_fun(mod_lat, slope_lat) %>% # FANCY FUNCTION
  mutate(dummy = "Absolute latitude")

# Predictions for scatter plot
# Generate values
lat_values <- new_df(abs_lat, meta_df) %>% # Function to create df with dummy vars
  rename(abs_lat = var) %>%
  mutate("abs_lat:ppCoral" = abs_lat*ppCoral,
         "abs_lat:ppEndosymbiont" = abs_lat*ppEndosymbiont,
         "abs_lat:ppPhytoplankton" = abs_lat*ppPhytoplankton,
         "abs_lat:ppSeagrass" = abs_lat*ppSeagrass) %>%
        # added_N_center = 0) %>%
  as.matrix()

# new model with non standardized background nutrients
model_lat_untran <- rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = meta_df,
                           method = "REML",
                           mods = ~ abs_lat*pp, # + added_N_center,
                           random = ~ 1 | paper_num/indiv, 
                           test = "t")

# use the generated values to plot background nutrients effect for each pp_group, nat_group held to 1
lat_predict <- as.data.frame(predict(object = model_lat_untran,
                                     newmods = lat_values)) %>%
  cbind(lat_values) %>%
  mutate(pp = c(rep("Algae", each = 100), 
                rep("Coral", each = 100), 
                rep("Endosymbiont", each = 100), 
                rep("Phytoplankton", each = 100), 
                rep("Seagrass", each = 100)))

# slope plot
slope_l <- int_effect_plot(int_effect_sum_l) +
  scale_colour_manual(values = rev(pal5)) +
  xlab("Slope")

# scatter plot
scatter_l <- scatter_fun(lat_predict, meta_df, abs_lat, rev(pal5)) + 
  labs(x = "Absolute latitude") +
  theme(axis.title.y = element_text(vjust = -15))

### Continuous background nitrogen --------------------------------------------

# Background nutrients
slope_bn <- coef(summary(mod_bn))[2,1] # extract the slope

# Calculate slopes for each pp
int_effect_sum <- slope_coeff_fun(mod_bn, slope_bn) %>% # FANCY FUNCTION
  mutate(dummy = "Background N")

# Predictions for scatter plot

# Generate the values
nut_values <- new_df(background_nut_umolL, meta_df_n) %>% # Function to create df with dummy vars
  rename(background_nut_umolL = var) %>%
  mutate("background_nut_umolL:ppCoral" = background_nut_umolL*ppCoral,
         "background_nut_umolL:ppEndosymbiont" = background_nut_umolL*ppEndosymbiont,
         "background_nut_umolL:ppPhytoplankton" = background_nut_umolL*ppPhytoplankton,
         "background_nut_umolL:ppSeagrass" = background_nut_umolL*ppSeagrass) %>%
       #  added_N_center = 0) %>%
  as.matrix()


# new model with non standardized background nutrients
model_back_untran <- rma.mv(yi = yi,
                            V = vi,
                            slab = author,
                            data = meta_df_n,
                            method = "REML",
                            mods = ~ background_nut_umolL*pp, # + added_N_center,
                            random = ~ 1 | paper_num/indiv, 
                            test = "t")

# use the generated values to plot background nutrients effect for each pp_group, nat_group held to 1
nut_predict <- as.data.frame(predict(object = model_back_untran,
                                     newmods = nut_values)) %>%
  cbind(nut_values) %>%
  mutate(pp = c(rep("Algae", each = 100), 
                rep("Coral", each = 100), 
                rep("Endosymbiont", each = 100), 
                rep("Phytoplankton", each = 100), 
                rep("Seagrass", each = 100)))

#seed
set.seed(123)

# slope plot
slope_n <- int_effect_plot(int_effect_sum) +
  scale_colour_manual(values = rev(pal5)) +
  xlab("Slope") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())

# scatter plot
scatter_n <- scatter_fun(nut_predict, meta_df_n,
                         background_nut_umolL, rev(pal5)) + 
  labs(x = expression("Background N"~(mu*M))) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())


### Continuous SST  -------------------------------------------
slope_sst <- coef(summary(mod_sst))[2,1] # extract the slope

# Calculate slopes
int_effect_sum_sst <- slope_coeff_fun(mod_sst, slope_sst) %>% # FANCY FUNCTION
  mutate(dummy = "Annual SST")

# Predictions for scatter plot
# Generate values
sst_values <- new_df(sst, meta_df) %>% # Function to create df with dummy vars
  rename(sst = var) %>%
  mutate("sst:ppCoral" = sst*ppCoral,
         "sst:ppEndosymbiont" = sst*ppEndosymbiont,
         "sst:ppPhytoplankton" = sst*ppPhytoplankton,
         "sst:ppSeagrass" = sst*ppSeagrass) %>%
       #  added_N_center = 0) %>%
  as.matrix()

# new model with non standardized background nutrients
model_sst_untran <- rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = meta_df,
                           method = "REML",
                           mods = ~ sst*pp,# + added_N_center,
                           random = ~ 1 | paper_num/indiv, 
                           test = "t")

# use the generated values to plot background nutrients effect for each pp_group, nat_group held to 1
sst_predict <- as.data.frame(predict(object = model_sst_untran,
                                     newmods = sst_values)) %>%
  cbind(sst_values) %>%
  mutate(pp = c(rep("Algae", each = 100), 
                rep("Coral", each = 100), 
                rep("Endosymbiont", each = 100), 
                rep("Phytoplankton", each = 100), 
                rep("Seagrass", each = 100)))

# slope plot
slope_sst <- int_effect_plot(int_effect_sum_sst, pal = rev(pal5)) +
  xlab("Slope") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  xlim(-1.77, 2.2)

# scatter plot
scatter_sst <- scatter_fun(sst_predict, meta_df_a, sst, rev(pal5)) + 
  labs(x = expression("Annual SST " ( degree*C))) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  xlim(-0.77, 30)

### Added N --------------------------------------------------------------------
slope_add <- coef(summary(mod_add))[2,1] # extract the slope

# Calculate slope for each pp
int_effect_sum_a <- slope_coeff_fun(mod_add, slope_add) # FANCY FUNCTION

# Predictions for scatter plot
# Generate values
add_values <- new_df(added_N_log, meta_added) %>% # Function to create df with dummy vars
  rename(added_N_log = var) %>%
  mutate("added_N_log:ppCoral" = added_N_log*ppCoral,
         "added_N_log:ppEndosymbiont" = added_N_log*ppEndosymbiont,
         "added_N_log:ppPhytoplankton" = added_N_log*ppPhytoplankton,
         "added_N_log:ppSeagrass" = added_N_log*ppSeagrass) %>%
  as.matrix()

# use the generated values to plot background nutrients effect for each pp_group, nat_group held to 1
add_predict <- as.data.frame(predict(object = mod_add,
                                     newmods = add_values)) %>%
  cbind(add_values) %>%
  mutate(pp = c(rep("Algae", each = 100), 
                rep("Coral", each = 100), 
                rep("Endosymbiont", each = 100), 
                rep("Phytoplankton", each = 100), 
                rep("Seagrass", each = 100)))

# slope plot
slope_a <- int_effect_plot(int_effect_sum_a) +
  scale_colour_manual(values = rev(pal5)) +
  xlab("Slope") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())

# scatter plot
scatter_a <- scatter_fun(add_predict, meta_added, added_N_log, rev(pal5)) + 
  labs(x = "Log(added N)") +
  theme(axis.title.y = element_text(vjust = 0))+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  xlim(-4.55, 4)


### Fig.4 plot 
# lat, sst, bn, added_N
(slope_l + slope_sst + slope_n +slope_a +
    plot_layout(ncol = 4)) / (scatter_l + scatter_sst + scatter_n + scatter_a +
                                plot_layout(ncol = 4)) +
  plot_annotation(tag_levels = "a")

# ggsave("Figures/Fig.4.cont_mods.png", device = "png", height = 125, width = 200, units = c("mm"), dpi = 600)


## Fig.5 Cat covariate plots --------------------------------------------------------------
# use the huge df we made for fig.3 again!

# plot the huge matrix of moderators x pp groups
ggplot(mega_df %>% filter(type != "All"), aes(emmean, variable, size = point_size)) +
  geom_point() +
  scale_size_manual(values = c("small" = 1, "normal" = 1.3)) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  facet_grid(
    rows = vars(type),
    cols = vars(pp),
    scales = "free_y",
    space = "free_y") +
  pub_theme() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)") +
  #scale_color_manual(values = (pal)) +
  # Add text labels inside the plots
  geom_text(aes(x = text_pos, label = n), # change to 1 if i add tt
            color = "black", size = 3, vjust = -0.5, hjust = 0)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

# ggsave("Figures/Fig.5.covariates.png", device = "png", height = 170, width = 210, units = c("mm"), dpi = 600)



## Fig.6 - Duration plots --------------------------------------------------
slope_dur <- coef(summary(mod_dur))[2,1] # extract the slope

# Calculate slope for each pp
int_effect_sum_d <- slope_coeff_fun(mod_dur, slope_dur) # FANCY FUNCTION

# Predictions for scatter plot
# Generate values
dur_values <- new_df(duration_log, meta_df %>% drop_na(duration_log)) %>% # Function to create df with dummy vars
  rename(duration_log = var) %>%
  mutate("duration_log:ppCoral" = duration_log*ppCoral,
         "duration_log:ppEndosymbiont" = duration_log*ppEndosymbiont,
         "duration_log:ppPhytoplankton" = duration_log*ppPhytoplankton,
         "duration_log:ppSeagrass" = duration_log*ppSeagrass) %>%
  as.matrix()

# use the generated values to plot background nutrients effect for each pp_group
dur_predict <- as.data.frame(predict(object = mod_dur,
                                     newmods = dur_values)) %>%
  cbind(dur_values) %>%
  mutate(pp = c(rep("Algae", each = 100), 
                rep("Coral", each = 100), 
                rep("Endosymbiont", each = 100), 
                rep("Phytoplankton", each = 100), 
                rep("Seagrass", each = 100)))

# slope plot
slope_d <- int_effect_plot(int_effect_sum_d) +
  scale_colour_manual(values = rev(pal5)) +
  xlab("Slope")

# scatter plot
scatter_d <- scatter_fun(dur_predict, meta_df, duration_log, rev(pal5)) + 
  labs(x = "Log(duration days)") +
  theme(axis.title.y = element_text(vjust = -15))

# Put it all together
(slope_d / scatter_d +
    plot_layout(ncol = 1)) +
  plot_annotation(tag_levels = "a")
# ggsave("Figures/Fig.6.dur.png", device = "png", height = 125, width = 90, units = c("mm"), dpi = 600)


# Supplemental Figures ---------------------------------------------------------

## Fig. S2 - LME map --------------------------------------------------
# Load the large marine ecosystems map
lme_m <- st_read("Data/lme_shp/lme66.shp") %>%
  left_join(meta_df %>%
              rename(LME_NAME = lme_name)%>%
              count(LME_NAME), by = "LME_NAME")

# make centroids 
lme_c <- lme_m %>%
  st_centroid(of_largest_polygon = TRUE) %>% # centroid becomes geometry
  st_transform(crs = "+proj=moll") # project to Mollweide

# Clip to a reasonable world extent (avoid actual corners)
clip <- st_as_sfc(st_bbox(c(xmin=-180, xmax=180, ymin=-85, ymax=85), crs=4326))

# clip the lme polys to that area
poly_clipped <- st_intersection(lme_m, clip)
poly_moll <- st_transform(poly_clipped, "+proj=moll")

# Clip world and points to the new elliptical mask
lme_clipped <- st_intersection(poly_moll, globe_mask_sf)

# set seed
set.seed(987654321)

# Plot with Mollweide and correct ellipse proportions
ggplot() +
  geom_sf(data = globe_mask_sf, fill = blue, color = NA) +
  geom_sf(data = world_clipped, fill = "white", color = "black") +
  # LME polygons
   geom_sf(data = lme_clipped, aes(fill = LME_NAME), color = "black", alpha = 0.5) +
  # add coordinates for each study
  geom_sf(data = points_clipped, pch = 21, colour = "black", fill = pal1, size = 2.5, alpha = 0.6) +
  # Add the black circle (ellipse) around the Earth
  geom_sf(data = globe_mask_ellipse, fill = NA, color = "black", size = 0.5) +  # Black circle outline
  # Circular labels for LMEs
    geom_sf_label(
      data = lme_c,
      aes(label = LME_NUMBER),   # replace with your label column
      fill = "white",
      color = "black",
      label.r = unit(0.5, "lines"),  # rounded corners
      label.size = 0.1,
      size = 1.1) +
  # Set the coordinate system to Mollweide and let it manage aspect ratio
  coord_sf(crs = st_crs("+proj=moll"), 
           expand = FALSE) +  # No manual scaling, just preserve aspect ratio
  pub_theme() +
  theme(
    axis.text = element_blank(),    # Remove axis labels
    axis.title = element_blank(),   # Remove axis titles
    axis.ticks = element_blank(),    # Remove axis ticks
    axis.line = element_blank(),    # Remove the axis lines
    panel.border = element_blank(),  # Remove the black border around the plot
    legend.position = "none"
  ) +
  scale_fill_manual(values = viridis::viridis(66))
# ggsave("Figures/Fig.S2_lme_map.png", width = 7.08661, height = 3, dpi = 600)


## Fig.S3 - LME effect sizes --------------------------------------------------
lme_m <- st_read("Data/lme_shp/lme66.shp") %>%
  left_join(
    meta_df %>%
      rename(LME_NAME = lme_name) %>%
      count(LME_NAME),
    by = "LME_NAME") %>%
  mutate(
    centroid = st_centroid(st_transform(geometry, 3857)) %>% 
      st_transform(4326),
    centroid_y = st_coordinates(centroid)[, 2]) %>%
  st_drop_geometry() %>%
  transmute(
    lme_name = LME_NAME,
    centroid_y)

# calculate centroid for "offshore"
offshore <- meta_df %>%
  filter(lme_name == "Offshore") %>%
  summarise(centroid_y = mean(lat_dec))

# make table
pred_lme_table <- as.data.frame(emmeans(prep_lme, specs="lme_name", type="response", tran = NA)) %>%
  left_join(meta_df %>%
              group_by(lme_name) %>%
              summarise(n = n()), by = "lme_name") %>%
  ungroup() %>%
  left_join(lme_m,
            by = "lme_name") %>%
  mutate(centroid_y = ifelse(lme_name == "Offshore", offshore$centroid_y, centroid_y),
         text_pos = case_when(
           lme_name == "Scotian Shelf" ~ emmean - 0.1,
           TRUE ~ as.numeric(emmean)),
         abs_lat = abs(centroid_y), # creating an absolute value lat
         tt = as.factor(ifelse(abs_lat >= 30, "Temperate", "Tropical")),
         climate = case_when(abs_lat <= 23.5 ~ "Tropical",
                             abs_lat <= 35 ~ "Subtropical",
                             TRUE ~ as.character("Temperate")))

# plot
ggplot(pred_lme_table, aes(emmean, reorder(lme_name, centroid_y))) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  pub_theme() +
  # theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)") +
  # Add text labels inside the plots
  geom_text(aes(x = text_pos, label = n), # change to 1 if i add tt
            color = "black", size = 3, vjust = -0.5, hjust = 0)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

# ggsave("Figures/Fig.S3.lme.png", device = "png", height = 160, width = 168, units = c("mm"), dpi = 600)


## Fig.S4 - Cat mods ----------------------------------------------------------
# make mega df
mega_df2 <- pred_bn_table1 %>%
  mutate(type = "N level") %>%
  rename(variable = bn_cat) %>%
  # add tt
  rbind(pred_tt_table1 %>%
          mutate(type = "Climate") %>%
          rename(variable = tt)) %>%
  # add sst
  rbind(pred_sst_table %>%
          mutate(type = "SST") %>%
          rename(variable = sst_cat)) %>%
  # tidy variables
  mutate(variable = case_when(variable == "Low nut" ~ "Low N",
                              variable == "High nut" ~ "High N",
                              TRUE ~ as.character(variable)),
         # nudge text position
         text_pos = case_when(
           pp == "Seagrass" & 
             type == "N level" & 
             variable == "High N" ~ emmean + 0.1,
           pp == "Seagrass" & 
             type == "SST" & 
             variable == "Cold" ~ emmean + 0.1,
           TRUE ~ emmean))

# plot
ggplot(mega_df2, aes(emmean, variable)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  facet_grid(
    rows = vars(type),
    cols = vars(pp),
    scales = "free_y",
    space = "free_y") +
  pub_theme() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)") +
  #scale_color_manual(values = (pal)) +
  # Add text labels inside the plots
  geom_text(aes(x = text_pos, label = n), # change to 1 if i add tt
            color = "black", size = 3, vjust = -0.5, hjust = 0) +
  #xlim(c(-1, 1.5)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3),
                     limits = c(-0.9,1.5))

# ggsave("Figures/Fig.S4.cat.png", device = "png", height = 70, width = 168, units = c("mm"), dpi = 600)


## Fig.S5 - Funnel plots ------------------------------------------------------
# For models involving moderators, the plot shows the residuals on the x-axis against their corresponding standard errors.
# residuals = ‘observed - fitted’ values

funnel.rma(mod_lat,
           xlab = "Residuals (log-scale)",
           atransf=exp,
           digits=2L,
           ylim = c(0, 1.24))

funnel.rma(mod_bn,
           xlab = "Residuals (log-scale)",
           atransf=exp,
           digits=2L,
           ylim = c(0, 1.24))

funnel.rma(mod_sst,
           xlab = "Residuals (log-scale)",
           atransf=exp,
           digits=2L,
           ylim = c(0, 1.24))
# A bit of a right skew

funnel.rma(mod_tt,
           xlab = "Residuals (log-scale)",
           atransf=exp,
           digits=2L,
           ylim = c(0, 1.24))

funnel.rma(mod_bn_cat,
           xlab = "Residuals (log-scale)",
           atransf=exp,
           digits=2L,
           ylim = c(0, 1.24))

funnel.rma(mod_sst_cat,
           xlab = "Residuals (log-scale)",
           atransf=exp,
           digits=2L,
           ylim = c(0, 1.24))
        #   at = c(0.05, 0.37, 2.72, 20.09))

funnel.rma(mod_add,
           xlab = "Residuals (log-scale)",
           atransf=exp,
           digits=2L,
           ylim = c(0, 1.24))

funnel.rma(mod_po4_pp,
           xlab = "Residuals (log-scale)",
           atransf=exp,
           digits=2L,
           ylim = c(0, 1.24))


## Fig.S7 - N supply vs c group -------------------------------------------------
sum_added_c <- ggpredict(mod_lm_add_c, terms = "c_group") %>% 
  rename(c_group = x,
         added_N = predicted)

# make palette that's smooth and nice (this is unnecessarily complicated)
# try again
library(scales)

# Reverse mako and combine
div_pal <- c(rev(viridisLite::viridis(100, option = "rocket")[14:100]), 
             viridisLite::viridis(100, option = "mako")[14:100])

# what % of the total range does the + and - sides of the spectrum make up???
max <- max(meta_added$yi) # 3.877808
min <- min(meta_added$yi) # -1.12493
range <- max - min # 5.002738

# calculate what % of the total spread each end makes
lower_range <- abs(min)/range # 0.2248628
upper_range <- abs(max)/range # 0.7751372

# so i want to find out what 22% of the lower half is?
(length(div_pal)/2)*lower_range # 19.56
# so that's the amount of units below 87 i want to go
(length(div_pal)/2) - 20 # so i want to start this thing at 67 i think????

# make pal
div_pal2 <- div_pal[67:172]

# make plot!
ggplot() +
  geom_jitter(data = meta_added, 
              aes(x = c_group, y = added_N, 
                  colour = yi, 
                  fill = yi,
                  pch = pp), 
              size = 1.5, alpha = 0.6, height = 0.2, width = 0.2)+
  # add dot + whiskers
  geom_point(data = sum_added_c,
             aes(x = c_group, y = added_N),
             size = 2) +
  geom_errorbar(data = sum_added_c,
                aes(x = c_group, 
                    y = added_N,
                    ymin = conf.low,
                    ymax = conf.high),
                width = 0.1,
                linewidth = 0.4)  +
  # clean up
  pub_theme()  +
  scale_shape_manual(values = 21:25) +
  scale_colour_gradientn(
    colours = div_pal2) +
  scale_fill_gradientn(
    colours = div_pal2) +
  labs(y = expression(paste("N supply rate"~(mu*M))),
       x = "Consumer group",
       colour = "LnRR",
       fill = "LnRR",
       pch = "Primary producer") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
# ggsave("Figures/Fig.S7.png", device = "png", height = 100, width = 140, units = c("mm"), dpi = 600)


## Fig.S8 - Seagrass vs epiphytes cont ----------------------------------------------
pal12 <- viridis::viridis(12)
pal_sge <- c(pal12[2], pal12[7])

slope_lat_sge <- coef(summary(mod_lat_sge))[2,1] # extract the slope

# Calculate slope for each pp
int_effect_sum_l_sge <- slope_coeff_fun_sge(mod_lat_sge, slope_lat_sge)

# Predictions for scatter plot
lat_values_sge <- new_df_sge(abs_lat, meta_sge) %>% # Function to create df with dummy vars
  rename(abs_lat = var) %>%
  mutate("abs_lat:ppSeagrass" = abs_lat*ppSeagrass) %>%
  as.matrix()

# new model with non standardized background nutrients
model_lat_untran_sge <- rma.mv(yi = yi,
                               V = vi,
                               slab = author,
                               data = meta_sge,
                               method = "REML",
                               mods = ~ abs_lat + pp + abs_lat:pp,
                               random = ~ 1 | paper_num/indiv, 
                               test = "t")

# use the generated values to plot background nutrients effect for each pp_group, nat_group held to 1
lat_predict_sge <- as.data.frame(predict(object = model_lat_untran_sge,
                                         newmods = lat_values_sge)) %>%
  cbind(lat_values_sge) %>%
  mutate(pp = c(rep("Algae", each = 100), 
                rep("Seagrass", each = 100)))

# slope plot
slope_l_sge <- int_effect_plot_sge(int_effect_sum_l_sge, pal_sge) +
  xlab("Slope")

# scatter plot
scatter_l_sge <- scatter_fun(lat_predict_sge, meta_sge, abs_lat, rev(pal_sge)) + 
  labs(x = "Absolute latitude") +
  theme(axis.title.y = element_text(vjust = -15))


# plots for sst
slope_sst_sge <- coef(summary(mod_sst_sge))[2,1] # extract the slope

# Calculate slope for each pp
int_effect_sum_sst_sge <- slope_coeff_fun_sge(mod_sst_sge, slope_sst_sge)

# Predictions for scatter plot
sst_values_sge <- new_df_sge(sst, meta_sge) %>% # Function to create df with dummy vars
  rename(sst = var) %>%
  mutate("sst:ppSeagrass" = sst*ppSeagrass) %>%
  as.matrix()

# new model with non standardized background nutrients
model_sst_untran_sge <- rma.mv(yi = yi,
                               V = vi,
                               slab = author,
                               data = meta_sge,
                               method = "REML",
                               mods = ~ sst + pp + sst:pp,
                               random = ~ 1 | paper_num/indiv, 
                               test = "t")

# use the generated values to plot background nutrients effect for each pp_group, nat_group held to 1
sst_predict_sge <- as.data.frame(predict(object = model_sst_untran_sge,
                                         newmods = sst_values_sge)) %>%
  cbind(sst_values_sge) %>%
  mutate(pp = c(rep("Algae", each = 100), 
                rep("Seagrass", each = 100)))

# slope plot
slope_sst_sge <- int_effect_plot_sge(int_effect_sum_sst_sge, pal_sge) +
  xlab("Slope")+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())

# scatter plot
scatter_sst_sge <- scatter_fun(sst_predict_sge, meta_sge, sst, rev(pal_sge)) + 
  labs(x = "SST") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())


# SAVE CONT MODS SEPERATELY
(slope_l_sge + slope_sst_sge +
    plot_layout(ncol = 2)) / (scatter_l_sge + scatter_sst_sge +
                                plot_layout(ncol = 2)) +
  plot_annotation(tag_levels = "a")

# ggsave("Figures/Fig.S8_sge_cont.png", device = "png", height = 125, width = 110, units = c("mm"), dpi = 600)


## Fig.S9 - seagrass vs epi cat ------------------------------------------------
# make predictions for bn
pred_table_sge_bn <- as.data.frame(emmeans(prep_sge_bn, specs="bn_cat", by = "pp", type="response", tran = NA)) %>% # emmeans to pull estimates
  mutate(y_num = 1:length(df)) %>%
  left_join((meta_sge %>%
               group_by(pp, bn_cat) %>%
               summarise(n = n())), by = c("pp", "bn_cat")) %>%
  ungroup()  %>%
  mutate(type = "N level") %>%
  rename(variable = bn_cat) %>%
  mutate(variable = case_when(variable == "Low nut" ~ "Low N",
                              variable == "High nut" ~ "High N",
                              TRUE ~ as.character(variable)),
         # nudge text position
         text_pos = case_when(
           pp == "Seagrass" & 
             variable == "High N" ~ emmean - 0.13,
           TRUE ~ emmean))

# plot
sge_bn_cat_plot <- 
  ggplot(pred_table_sge_bn, aes(emmean, variable)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  facet_grid(
    cols = vars(pp),
    scales = "free_y",
    space = "free_y") +
  pub_theme() +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)") +
  # Add text labels inside the plots
  geom_text(aes(x = text_pos, label = n), # change to 1 if i add tt
            color = "black", size = 3, vjust = -0.5, hjust = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3),
                     limits = c(-1, 1.5))

sge_bn_cat_plot
# ggsave("Figures/Fig.S9_sge_cat.png", device = "png", height = 70, width = 80, units = c("mm"), dpi = 600)


## Fig.S2.1 - 2.2 Continuous background phosphate ------------------------------------

# Background phosphate
slope_po4 <- coef(summary(mod_po4_pp))[2,1] # extract the slope

# Calculate slopes for each pp
int_effect_sum <- coef(summary(mod_po4_pp)) %>%
  slice(-(1)) %>%
  slice(-(2:5)) %>%
  mutate(int = factor(c("Algae", "Coral", "Phytoplankton", "Seagrass")),
         effect = ifelse(int != "Algae", estimate + slope_po4,
                         slope_po4), # add intercept to the estimated coefficient
         ci.lb_adjust = effect - (1.96*se), #adjust the CI based on the new estimate
         ci.up_adjust = effect + (1.96*se),
         dummy = 1:4) %>% 
  mutate(dummy = "Background P") %>%
  bind_rows(tibble(int = "Endosymbiont"))


# Predictions for scatter plot
# Make new DF to generate predictions from
meta_df_p2 <- meta_df_p %>% 
  filter(pp != "Endosymbiont") %>%
  drop_na(po4)

# Generate the values
mm <- meta_df_p2 %>%
  group_by(pp) %>%
  summarise(min = min(po4),
            max = max(po4)) %>%
  mutate(x = c(1,2,3,4)) %>%
  as.data.frame()

nut_values <- data.frame(
  var = c(
    seq(from = mm[1,2], to = mm[1,3], length.out = 100),
    seq(from = mm[2,2], to = mm[2,3], length.out = 100),
    seq(from = mm[3,2], to = mm[3,3], length.out = 100),
    seq(from = mm[4,2], to = mm[4,3], length.out = 100)),
  ppCoral = c(rep(0, each = 100), 
              rep(1, each = 100), 
              rep(0, each = 100), 
              rep(0, each = 100)),
  ppPhytoplankton = c(rep(0, each = 100), 
                      rep(0, each = 100),
                      rep(1, each = 100), 
                      rep(0, each = 100)),
  ppSeagrass = c(rep(0, each = 100), 
                 rep(0, each = 100), 
                 rep(0, each = 100), 
                 rep(1, each = 100))) %>%
  rename(po4 = var) %>%
  mutate("po4:ppCoral" = po4*ppCoral,
         #   "po4:ppEndosymbiont" = po4*ppEndosymbiont,
         "po4:ppPhytoplankton" = po4*ppPhytoplankton,
         "po4:ppSeagrass" = po4*ppSeagrass) %>%
  as.matrix()


# new model with non standardized background nutrients
model_po4_untran <- rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = meta_df_p2,
                           method = "REML",
                           mods = ~ po4*pp,
                           random = ~ 1 | paper_num/indiv, 
                           test = "t")

# use the generated values to plot Background Putrients effect for each pp_group, nat_group held to 1
po4_predict <- as.data.frame(predict(object = model_po4_untran,
                                     newmods = nut_values)) %>%
  cbind(nut_values) %>%
  mutate(pp = c(rep("Algae", each = 100), 
                rep("Coral", each = 100), 
                #  rep("Endosymbiont", each = 100), 
                rep("Phytoplankton", each = 100), 
                rep("Seagrass", each = 100)))

#seed
set.seed(123)

# slope plot
slope_p <- int_effect_plot(int_effect_sum) +
  scale_colour_manual(values = rev(pal5)) +
  xlab("Slope") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())

# scatter plot
scatter_p <- scatter_fun(po4_predict, meta_df_p,
                         po4, rev(pal5)) + 
  labs(x = expression("Background P"~(mu*M))) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "right")
slope_p + scatter_p

# plot predictions facet wrapped
ggplot() +
  geom_point(data = meta_df_p, 
             aes(po4, yi, colour = pp, fill = pp), 
             size = 1.2, alpha = 0.5) +
  geom_line(data = po4_predict, aes(po4, pred, colour = pp),  linewidth = 0.5) +
  geom_ribbon(data = po4_predict, aes(po4, pred, ymin = ci.lb, ymax = ci.ub, fill = pp), 
              alpha = 0.1) +
  labs(y = "Effect (LnRR)") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.2) +
  scale_colour_manual(values = rev(pal5)) +
  scale_fill_manual(values = rev(pal5)) +
  pub_theme() +
  theme(legend.position = "none") +
  labs(colour = "", fill = "")  +
  facet_wrap(~pp)+ 
  labs(x = expression("Background P"~(mu*M)))

# ggsave("Figures/Fig.S2.1.po4.png", device = "png", height = 125, width = 200, units = c("mm"), dpi = 600)

# just phyto
ggplot(data = meta_df_p %>% filter(pp == "Phytoplankton"), 
       aes(po4, yi)) +
  geom_point(size = 1.2, alpha = 0.5, colour ="#2A788EFF") +
  geom_smooth(method = "lm", colour ="#2A788EFF", fill ="#2A788EFF", alpha = 0.3) +
  labs(y = "Effect (LnRR)") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.2) +
  pub_theme() +
  theme(legend.position = "none") +
  labs(x = expression("Background P"~(mu*M)))

# ggsave("Figures/Fig.S2.2.po4.phyto.png", device = "png", height = 90, width = 90, units = c("mm"), dpi = 600)


# Model selection --------------------------------------------------------------
mod_list <- c("lat_stand", "tt",
              "bn_stand", "bn_cat",
              "sst_stand", "sst_cat")

covar_list <- c("nat_bin2", "grow_time_end", "res", "c_group", 
                "duration_stand", "added_N_stand")

combos <- crossing(moderator = mod_list, covar = covar_list)

results <- apply(combos, 1, function(row) {
  meta_aic(row["moderator"], row["covar"], meta_df)
}, simplify = FALSE)

# bind results into df!!!
results_df <- bind_rows(results)

# Replacement mapping
replacements <- c(
  "lat_stand" = "lat",
  "tt" = "lat cat",
  "bn_stand" = "BN",
  "bn_cat" = "BN cat",
  "sst_stand" = "SST",
  "sst_cat" = "SST cat",
  "nat_bin2" = "ED",
  "grow_time_end" = "measure",
  "res" = "response",
  "c_group" = "C taxa",
  "duration_stand" = "duration",
  "added_N_stand" = "added N"
)

# Apply replacements in tidyverse style
results_df_clean <- results_df %>%
  mutate(model = str_replace_all(model, replacements))

# write_csv(results_df_clean, "Tables/covar_aic_results2.csv")


# More selection
aic <- AIC(mod_null, mod_lat, mod_tt, mod_bn, mod_sst, mod_sst_cat, mod_nat, mod_nat_pp, mod_inst_grow, mod_inst_grow_pp, mod_res, mod_res_pp, mod_taxa, mod_taxa_pp, mod_pp, mod_dur, mod_add, mod_lme)


# In general, all three ways I’ve tried to do these stats give a similar answer.

# Correlation between variables
cor(cbind(meta_data$lat_stand, meta_data$bn_stand, meta_data$temp_stand), use = "complete.obs")
# so temp, lat, and background_nutrients are correlated (|value| < 0.5)


# Summary stats ----------------------------------------------------------------

# climate regions
meta_data %>% count(tt)

# How many pp were there in each group?
meta_data %>% 
  count(pp)

# table of pp species used
pp_used <- meta_data %>%
  count(pp, pp_family, pp_genus, pp_species) # typo in datasheet
pp_used
# write_csv(pp_used, "Tables/pp_used.csv")

# Indiv, pop, performance responses
d <- meta_data %>%
  count(meas_type, response)

# primary producer responses measured
metrics <- meta_data %>%
  count(response) %>%
  arrange(-n)

# But what did they actually measure for each response
metrics2 <- meta_data %>%
  count(pp, response) %>%
  rename(small_n = n) %>%
  left_join(metrics, by = "response")
# write_csv(metrics, "Tables/metrics.csv")

# how many experiments were in each experimental design category?
nat_counts <- meta_data %>% 
  count(nat_bin2) %>%
  mutate(percent = 100*n/length(meta_data$obs_num))

# how many things measured for each pp study
metrics2 <- meta_data %>%
  group_by(pp, paper_num, year_of_study) %>%
  count(response) %>%
  summarise(n=n())


# Consumers 
# % of consumers used in each broad category?
meta_data %>%
  count(c_activity)

# which taxa were in each group?
consumer_taxa <- meta_data %>%
  count(c_taxa, c_taxa)

# Table of consumer species used
c_used <- meta_data %>%
  count(c_group, c_family, c_genus, c_species)
# write_csv(c_used, "Tables/c_used.csv")
  
# Cross of consumer x pp species used
cross <- meta_data %>%
  count(pp, c_taxa) %>%
  mutate(percent = 100*n/sum(n))


## Proxy info 

# background nutrients
bn_proxy <- meta_data %>%
  as.data.frame() %>%
  filter(bn_proxy == "proxy") %>%
  transmute(author = author,
            pub_year = year,
            `Year of study` = year_of_study,
            `Proxy year` = bn_proxy_date,
            `Proxy distance (km)` = bn_proxy_dist) %>%
  unique() %>%
  arrange(-desc(author))
# write_csv(bn_proxy, "Tables/bn_proxy.csv")


temp_proxy <- meta_data %>%
  as.data.frame() %>%
  filter(temp_proxy == "proxy")  %>%
  transmute(author = author,
            pub_year = year,
            `Year of study` = year_of_study,
            `Proxy year` = temp_proxy_date,
            `Proxy distance (km)` = temp_proxy_dist_km) %>%
  unique() %>%
  arrange(-desc(author))
# write_csv(temp_proxy, "Tables/temp_proxy.csv")


# more classification
bn_classification <- meta_data %>%
  group_by(author, year, nutrient_levels, bn_cat, bn_proxy) %>%
  summarize(bn = mean(background_nut_umolL))%>%
  arrange(bn) %>%
  ungroup() %>%
  mutate(Paper = word(author, 1),
         Paper = trimws(Paper, whitespace = ","),
         Paper = sub("$", " et al.,", Paper), # $ specifies to add the text to the end
         Paper =  paste(Paper, year)) %>%
  transmute(Paper = Paper,
            author = author, # keep author so I can fix issues in excel
            `Author description` = nutrient_levels, 
            `Proxy or measured` = bn_proxy,
            `Nitrogen (uM)` = round(bn, 2),
            `Nitrogen level` = word(bn_cat, 1)
            ) 

# write_csv(bn_classification, "Tables/bn_classification.csv")

# Giant table of all the data

huge_data_table <- meta_data %>%
  mutate(Paper = word(author, 1),
         Paper = trimws(Paper, whitespace = ","),
         Paper = sub("$", " et al.,", Paper), # $ specifies to add the text to the end
         Paper =  paste(Paper, year)) %>%
  count(Paper, author, nat_bin2, c_group, pp, duration_days, grow_time_end, res, bn_cat, tt, added_N) %>%
  arrange(Paper) %>%
  transmute(author = author, 
            Paper = Paper,
            Climate = case_when(tt == "Temperate" ~ "temp",
                                tt == "Tropical" ~ "trop"),
            `N levels` = word(bn_cat, 1),
            Producer = case_when(pp == "Endosymbiont" ~ "Endo",
                                 pp == "Phytoplankton" ~ "Phyto",
                                 TRUE ~ as.character(pp)), 
            Design = case_when(nat_bin2 == "Unstocked"  ~ "unstocked",
                               nat_bin2 == "Stocked"  ~ "stocked",
                               nat_bin2 == "Mesocosm"  ~ "mesocosm"),
            Measured = case_when(grow_time_end == "growth_rate_data" ~ "growth",
                                 grow_time_end == "endpoint" ~ "end point",
                                 grow_time_end == "timepoint" ~ "inst"),
            Response = res,
            Consumer = ifelse(c_group == "Sea cucumbers", "Cukes", as.character(c_group)),
            Duration = round(duration_days, 2),
            `Added N` = round(added_N, 2),
           # Duration = ifelse(duration_days == 0, NA, round(duration_days, 2)),
           # `Added N` = ifelse(added_N == 0, NA, round(added_N, 2)),
            n = n)

# write_csv(huge_data_table, "Tables/huge_data_table_bigger2.csv")


# Are individual and population measurements relatively evenly distributed across ecosystem, background nutrients, etc
response_tt <- meta_data %>%
  count(meas_type, tt) %>%
  group_by(meas_type) %>%
  mutate(percent = n/sum(n))
# for each measurement type, half of the effect sizes are from tropics vs temperates

response_bn <- meta_data %>%
  count(meas_type, bn_cat) %>%
  group_by(meas_type) %>%
  mutate(percent = n/sum(n))
# most performance measurements were in low nutrient regions

# what did they measure when they measured growth?
growth_metrics <- meta_data %>%
  filter(growth_rate == "Y") %>%
  mutate(response = case_when(response == "mass" ~ "biomass",
                              TRUE ~ as.character(response))) %>%
  count(response)

# What metrics were used when growth wasn't measured
other_metrics <- meta_data %>%
  filter(growth_rate == "N") %>%
  mutate(response = case_when(response == "mass" ~ "biomass",
                              TRUE ~ as.character(response))) %>%
  count(response)

# Look at individual pp metrics to sea what was measured for each
# Run this for each metric
measured <- meta_data %>%
  filter(response == "width") %>%
  select(response, measured, growth_rate, pp_taxa, comments)


# how many effect sizes for bn vs po4?
nrow(meta_df %>% filter(bn_proxy == "measure") %>% count(title)) # 127 effect sizes, 19 studies
nrow(meta_df %>% drop_na(po4)%>% count(title)) # 77, 13 studies
nrow(meta_df %>% filter(temp_proxy == "measure")%>% count(title)) # 168 effect sizes, 23 studies

nrow(meta_df %>% drop_na(added_N)) # 132
nrow(meta_df %>% drop_na(added_po4)) # 75

n_bn_pp <- meta_df %>% filter(bn_proxy == "measure") %>% count(pp)
n_po4_pp <- meta_df %>% drop_na(po4)%>% count(pp)
n_po4_add_pp <- meta_df %>% drop_na(added_po4)


# Independence and Bias --------------------------------------------------------

## Rosenthal's Failsafe -----
fsn(x = yi,  vi = vi, data = meta_data, type="Rosenthal") #file drawer analysis (missing studies have an effect size = 0)
# Fail-safe N: 99,557
# That seems like quite a large number?????


## Check for influential studies -----------------------------------------------
# ?metafor::influence.rma.mv()
# dfbetas(model_full_back, progbar=FALSE, reestimate=TRUE, parallel="no", ncpus=1, cl=NULL, ...)
# hatvalues(model_full_back, type="diagonal")

# Change reestimate to FALSE and these fun faster

# Good help here
# https://wviechtb.github.io/metafor/reference/influence.rma.mv.html

# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/heterogeneity.html?q=cook#cooks-distance

# Calculate cook's distance for the lat model
cook_lat <- cooks.distance(mod_lat, cluster = meta_df$paper_num, progbar=FALSE, reestimate=FALSE, parallel="no", ncpus=1, cl=NULL) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(paper_num = 1, cook_lat = 2)

# cooks for bn
meta_bn <- meta_df %>% filter(bn_proxy == "measure")
cook_bn <- cooks.distance(mod_bn, cluster = meta_bn$paper_num, progbar=FALSE, reestimate=FALSE, parallel="no", ncpus=1, cl=NULL) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(paper_num = 1, cook_bn = 2)

#sst
cook_sst <- cooks.distance(mod_sst, cluster = meta_df$paper_num, progbar=FALSE, reestimate=FALSE, parallel="no", ncpus=1, cl=NULL) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(paper_num = 1, cook_sst = 2)

# cooks for phosphate
cook_p <- cooks.distance(mod_po4, cluster = meta_df_p$paper_num, progbar=FALSE, reestimate=FALSE, parallel="no", ncpus=1, cl=NULL) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(paper_num = 1, cook_p = 2)

# cat mods
cook_lat_c <- cooks.distance(mod_tt, cluster = meta_df$paper_num, progbar=FALSE, reestimate=FALSE, parallel="no", ncpus=1, cl=NULL) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(paper_num = 1, cook_lat_c = 2)

# Calculate the cook's distance for the nutrient cat model 
meta_bn_c <- meta_df %>% drop_na(bn_cat)
cook_bn_c <- cooks.distance(mod_bn_cat, cluster = meta_bn_c$paper_num, progbar=FALSE, reestimate=FALSE, parallel="no", ncpus=1, cl=NULL) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(paper_num = 1, cook_bn_c = 2)

# cat mods
cook_sst_c <- cooks.distance(mod_sst, cluster = meta_df$paper_num, progbar=FALSE, reestimate=FALSE, parallel="no", ncpus=1, cl=NULL) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(paper_num = 1, cook_sst_c = 2)

# add cook's distance for these two models to the meta data
# dmetar uses cook's distance > 0.45 as the cut-off for influental studies
cook_data <- cook_lat %>%
  left_join(cook_bn, by = "paper_num") %>%
  left_join(cook_sst, by = "paper_num") %>%
  left_join(cook_lat_c, by = "paper_num") %>%
  left_join(cook_bn_c, by = "paper_num") %>%
  left_join(cook_sst_c, by = "paper_num") %>%
  left_join(cook_p, by = "paper_num") %>%
 # left_join(meta_data , by = "paper_num") %>%
  mutate(influential = ifelse(cook_lat > 0.45 | 
                                cook_bn > 0.45 |
                                cook_sst > 0.45 |
                                cook_p > 0.45, "influential", "ok"),
         problem = ifelse(cook_lat > 1 |
                            cook_bn > 1 |
                            cook_sst > 1 |
                            cook_p > 1, "problem", "ok")) %>%
  mutate_if(is.numeric, round, digits = 2)

# Is cook's distance similar for each?
# We're looking for a perfect straight line
# Not exactly but not terrible?
ggplot(cook_data, aes(cook_lat, cook_bn, label = paper_num, colour = influential)) +
  geom_jitter() +
  geom_text(nudge_x = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0.45, lty = "dashed")+
  geom_vline(xintercept = 0.45, lty = "dashed")

# Drop papers with lots of obs
# make a df where each row is a paper, and num_obs is the number of observations that paper contributed
max_obs <- meta_df %>%
  group_by(paper_num) %>%  
  summarize(n_obs = n_distinct(obs_num))
hist(max_obs$n_obs)
sum(max_obs$n_obs) # 246
# 16 = Slattery cave sponge paper has more than 10% of the observations


# rerun latitude
# mods that are high for lat: 22, 16, 27, 20, 40, 23
# Remove those papers one at a time and rerun the stats! 
meta_data_drop <- meta_df %>%
  filter(paper_num != 27) 
# filter(!(paper_num %in% c(22, 16, 27, 20, 40, 23)))

mod_lat_redo <- rma.mv(yi = yi,
                      V = vi,
                      slab = author,
                      data = meta_data_drop,
                      method = "REML",
                      mods = ~ lat_center*pp,
                      random = ~ 1 | paper_num/indiv, 
                      test = "t")
summary(mod_lat_redo)
# 22 coral and endo int = signif, seagrass more diff
# 16, 40 coral:lat interaction becomes signif
# 20 lat slope goes away, no signif interactions
# 23 endo int becomes signif

# Original summary
# estimate      se     tval   df    pval    ci.lb    ci.ub      
# intrcpt                       0.4331  0.1262   3.4311  236  0.0007   0.1844   0.6818  *** 
# lat_center                    0.0211  0.0099   2.1284  236  0.0343   0.0016   0.0406    * 
# ppCoral                      -0.2021  0.3287  -0.6148  236  0.5393  -0.8496   0.4454      
# ppEndosymbiont               -0.5911  0.3718  -1.5898  236  0.1132  -1.3236   0.1414      
# ppPhytoplankton              -0.0157  0.1749  -0.0898  236  0.9285  -0.3602   0.3288      
# ppSeagrass                   -0.1957  0.1619  -1.2084  236  0.2281  -0.5146   0.1233      
# lat_center:ppCoral           -0.0383  0.0247  -1.5546  236  0.1214  -0.0869   0.0102      
# lat_center:ppEndosymbiont    -0.0472  0.0307  -1.5392  236  0.1251  -0.1076   0.0132      
# lat_center:ppPhytoplankton   -0.0286  0.0127  -2.2444  236  0.0257  -0.0537  -0.0035    * 
# lat_center:ppSeagrass        -0.0370  0.0121  -3.0620  236  0.0025  -0.0608  -0.0132   ** 


# rerun bn
# mods that are high for bn: 16, 42, 30, 19, 43, 18
meta_data_drop <- meta_df %>%
  filter(paper_num != 30)
#  filter(!(paper_num %in% c(16, 42, 30, 19, 43, 18)))

# rerun background nut mod
mod_bn_redo <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_data_drop %>% filter(bn_proxy == "measure"),
                          method = "REML",
                          mods = ~ bn_center*pp,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")
summary(mod_bn_redo)

# drop 16, bn goes away
# drop 30, bn_center:ppPhytoplankton becomes signif

# after moving to filtered df, dropping 30, 19, and 43 makes bn MORE signif

#                              estimate      se     tval   df    pval    ci.lb    ci.ub      
#  intrcpt                      0.6220  0.1785   3.4853  117  0.0007   0.2686   0.9755  *** 
#  bn_center                    0.0445  0.0222   2.0081  117  0.0469   0.0006   0.0884    * 
#  ppCoral                     -0.2074  0.9734  -0.2131  117  0.8316  -2.1352   1.7204      
#  ppEndosymbiont              -0.2065  0.7488  -0.2757  117  0.7832  -1.6894   1.2764      
#  ppPhytoplankton             -0.4527  0.2101  -2.1547  117  0.0332  -0.8687  -0.0366    * 
#  ppSeagrass                  -0.4872  0.2989  -1.6298  117  0.1058  -1.0793   0.1048      
#  bn_center:ppCoral           -0.1650  0.3086  -0.5348  117  0.5938  -0.7761   0.4461      
#  bn_center:ppEndosymbiont     0.0220  0.2340   0.0941  117  0.9252  -0.4414   0.4854      
#  bn_center:ppPhytoplankton   -0.0481  0.0355  -1.3537  117  0.1785  -0.1184   0.0223      
#  bn_center:ppSeagrass        -0.0383  0.1041  -0.3679  117  0.7136  -0.2446   0.1679      

# rerun sst
# mods that are high for sst: 23, 22, 16, 27, 20, 40
meta_data_drop <- meta_df %>%
  filter(paper_num != 19)
#  filter(!(paper_num %in% c(23, 22, 16, 27, 20, 40)))

# rerun sst mod
mod_sst_redo <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_data_drop,
                          method = "REML",
                          mods = ~ sst_center*pp,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")
summary(mod_sst_redo) # no diffs!!!!
# 22 seagrass becomes signif
#estimate      se     tval   df    pval    ci.lb    ci.ub      
#intrcpt                       0.4388  0.1189   3.6891  236  0.0003   0.2045   0.6731  *** 
#sst_center                   -0.0397  0.0153  -2.5896  236  0.0102  -0.0699  -0.0095    * 
#ppCoral                      -0.7558  0.6616  -1.1423  236  0.2545  -2.0593   0.5477      
#ppEndosymbiont               -0.0747  0.7545  -0.0991  236  0.9212  -1.5611   1.4117      
#ppPhytoplankton              -0.0163  0.1606  -0.1017  236  0.9191  -0.3328   0.3002      
#ppSeagrass                   -0.2326  0.1506  -1.5448  236  0.1237  -0.5293   0.0640      
#sst_center:ppCoral            0.1507  0.0943   1.5983  236  0.1113  -0.0351   0.3365      
#sst_center:ppEndosymbiont     0.0107  0.1098   0.0978  236  0.9221  -0.2056   0.2271      
#sst_center:ppPhytoplankton    0.0587  0.0201   2.9208  236  0.0038   0.0191   0.0984   ** 
#sst_center:ppSeagrass         0.0739  0.0206   3.5882  236  0.0004   0.0333   0.1144  *** 


# rerun po4
# mods that are high for po4: 33, 30, 26, 4, 11, 18, 9, and 42 is only endo paper
# BUT if we just include po4 (no pp) in model, its just 30 and 9
meta_data_drop <- meta_df %>%
  filter(paper_num != 9) #%>%
#  filter(!(paper_num %in% c(33, 26, 4, 18)))


# Model with continuous background phosphate
mod_po4_redo <- rma.mv(yi = yi,
                  V = vi,
                  slab = author,
                  data = meta_data_drop,
                  method = "REML",
                  mods = ~ po4_stand,
                  random = ~ 1 | paper_num/indiv, 
                  test = "t")
summary(mod_po4_redo)

# 33 intercept changes and becomes signif, po4 becomes negative slope, so many things
# 26 int becomes signif
# 4 and 18 coral and endo become signif

# simple model: without 30, po4 becomes non signif
# 

#                             estimate   se     tval    df    pval    ci.lb    ci.ub     
# intrcpt                      0.1041  0.2509   0.4151  68  0.6794  -0.3966   0.6049     
# po4_stand                    0.7704  0.2287   3.3682  68  0.0012   0.3140   1.2268  ** 
# ppCoral                      1.6036  0.9488   1.6901  68  0.0956  -0.2897   3.4970   . 
# ppEndosymbiont               0.5657  0.4359   1.2977  68  0.1988  -0.3042   1.4355     
# ppPhytoplankton              0.0406  0.2871   0.1414  68  0.8880  -0.5324   0.6136     
# ppSeagrass                   1.1526  0.8647   1.3329  68  0.1870  -0.5730   2.8782     
# po4_stand:ppCoral            1.6442  2.2532   0.7297  68  0.4681  -2.8520   6.1404     
# po4_stand:ppPhytoplankton   -0.6914  0.2406  -2.8740  68  0.0054  -1.1714  -0.2113  ** 
# po4_stand:ppSeagrass         1.3942  1.9201   0.7261  68  0.4703  -2.4373   5.2257     


# Exploratory stuff ------------------------------------------------------
# Mess around with excluding experiments that are too natural
# make a df without any fully just "there's a natural thing happening and we're just sampling it" papers
new <- meta_data %>%
  mutate(too_nat = case_when(grow_time_end == "endpoint" & 
                               nat_bin2 == "Unstocked" & 
                               response == "percent cover" ~ "too natural",
                             grow_time_end == "endpoint" & 
                               nat_bin2 == "Unstocked" & 
                               duration_comments == "just sampled natural eelgrass beds" ~ "too natural",
                             grow_time_end == "endpoint" & 
                               nat_bin2 == "Unstocked" & 
                               duration_comments == "just sampled natural eelgrass beds. no point in trying to convert this to after - before all metrics declined over the sampling period it would all be negative" ~ "too natural",
                             TRUE ~ as.character("ok")
  )) %>%
  filter(too_nat == "ok")


good <- meta_data %>%
  # want these
  filter(!(nat_bin == "Mesocosm" | 
             nat_bin == "Field: both stocked" |
             nat_bin == "Field: animals stocked" |
             nat_bin == "Field: PP stocked" |
             nat_bin == "Field: unstocked" & growth_rate == "Y")) %>%
  # don't want this
  filter(!(grow_time_end == "endpoint" & response == "percent cover")) 

meta_unstock <- meta_data %>%
  filter(nat_bin == "Field: unstocked") %>%
  mutate(level_manip = case_when(
    paper_num %in% c(3, 4, 8, 14, 37, 10, 11, 15, 19, 24, 32, 33) ~ "within site",
    paper_num %in% c(39, 41, 40, 16, 6, 22) ~ "amoung site"
  ),
  useable = case_when(level_manip == "within site" ~ "good",
                      level_manip == "amoung site" & growth_rate == "Y" ~ "good",
                      level_manip == "amoung site" & growth_rate == "N" ~ "bad")) %>%
  select(indiv, level_manip, useable)
# good ones (within site manipulation?) = paper num 3, 4, 8, 14, 37, 10, 11, 15, 19, 24, 32, 33,

# unsure ones
# unsure ones (seabird vs not but replicated, also measured growth) = 39, 41
# unsure ones (seabird vs not but replicated, but just measured % cover) = 40
# unsure ones (proximity to cave replicated, some % cover some growth) = 16

# healing, vinther, slattery, marmen are the papers that were amoung site, and also didn't measure growth

# bad ones
# 6 = seabird vs not (site level), just measured ambient bottom water chl a (no growth)
# 22 = mussel vs mixed site (site level, only one of each site type), measured ambient stuff and also growth

# what do we think about seabird islands??? i think those should count???? should they count when they measure growth? probably? but what about when they just measure ambient percent cover???

# I'm saying that artificial reef studies are good, and studies where they looked at corals with vs without fish are good
# cut out the marmen seabird bottom chl a sampling, the vinther paper where they just compared a non aquculture site to an aquaculture site, and that dumb slattery cave study???

# basically getting rid of observational studies where the entire site is the level of study???

meta_df3 <- meta_data %>% 
  left_join(meta_unstock, by = "indiv") %>%
  filter((useable != "bad") %>% replace_na(TRUE)) %>%
  filter(paper_num != "30")
# We lose Healing and Marmen completely, but Slattery and Vinther are retained in meta_df because they had growth measurements

# added n vs pp
w <- meta_df %>% group_by(c_taxa) %>% summarise(min = min)


## fish trophic levels ------
fish <- meta_df %>% filter(c_taxa == "F") %>%
  select(title, author, year, c_species_family_comm, c_family, c_genus, c_species) %>%
  unique() %>%
  mutate(trophic = case_when(
    #   c_genus == "Haemulon" ~ "invertivorous", # crustaceans?
    c_family == "Haemulidae" ~ "omnivore",
    c_genus == "Dicentrarchus" ~ "predator", # piscivorous
    c_genus == "Dascyllus" & c_species == "aruanus" ~ "omnivore", # zooplankton, benthic invertebrates, and algae
    c_genus == "Dascyllus" & c_species ==  "flavicaudus" ~ "invertivorous",
    c_genus == "Amphiprion" ~ "omnivore", # small invertebrates and algae
    c_genus == "Salmo" ~ "predator", # adults at sea feed on squids, shrimps, and fish
    c_genus == "Siganus" ~ "herbivore",
    c_species_family_comm == "community" ~ "community"
  ))

meta_fish <- meta_df %>% 
  filter(c_taxa == "F") %>%
  mutate(trophic = case_when(
 #   c_genus == "Haemulon" ~ "invertivorous", # crustaceans?
    c_family == "Haemulidae" ~ "omnivore",
    c_genus == "Dicentrarchus" ~ "predator", # piscivorous
    c_genus == "Dascyllus" & c_species == "aruanus" ~ "omnivore", # zooplankton, benthic invertebrates, and algae
    c_genus == "Dascyllus" & c_species ==  "flavicaudus" ~ "invertivorous",
    c_genus == "Amphiprion" ~ "omnivore", # small invertebrates and algae
    c_genus == "Salmo" ~ "predator", # adults at sea feed on squids, shrimps, and fish
    c_genus == "Siganus" ~ "herbivore",
    c_species_family_comm == "community" ~ "community"
    ))

mod_fish <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = meta_fish,
                 method = "REML",
                 mods = ~ trophic,
                 random = ~ 1 | paper_num/indiv, 
                 test = "t")
summary(mod_fish)

# prep model for emmeans
prep_fish <- emmeans_prep(mod_fish) # use function to prep model for emmeans

# pairwise
summary(pairs(emmeans(prep_fish,  ~ trophic)))

# make table
pred_fish_table <- as.data.frame(emmeans(prep_fish, specs="trophic", type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df)) %>%
  left_join((meta_fish %>%
               group_by(trophic) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("trophic")) %>% 
  ungroup()%>%
  mutate(pp = "All")

ggplot(pred_fish_table, aes(emmean, trophic)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  pub_theme() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)") +
  #scale_color_manual(values = (pal)) +
  # Add text labels inside the plots
  geom_text(aes(x = emmean, label = n), # change to 1 if i add tt
            color = "black", size = 3, vjust = -0.5, hjust = 0)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

# summary fish
meta_fish_summary <- meta_fish %>%
  select(title, trophic) %>%
  unique() %>%
  count(trophic)

# for fun
mod_bn <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = meta_df %>% filter(bn_proxy == "measure"),
                 method = "REML",
                 mods = ~ bn_center*tt,
                 random = ~ 1 | paper_num/indiv, 
                 test = "t")
summary(mod_bn)

ggplot(meta_df %>% filter(bn_proxy == "measure"), 
       aes(background_nut_umolL, yi)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = lm) +
  facet_wrap(~tt)

ggplot(meta_df, 
       aes(sst, yi)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = lm) +
  facet_wrap(~tt)

ggplot(meta_df, 
       aes(abs_lat, yi)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = lm) +
  facet_wrap(~tt)


# reviewer suggestion to look at sign of effect???
d <- meta_df %>%
  mutate(pos = ifelse(yi > 0, "positive", "negative")) %>%
  group_by(pp) %>%
    count(pos)
  
ggplot(d, aes(pos, n, colour = pos)) +
  geom_jitter() +
  facet_wrap(~pp)

# region by ocean?
mod_o <- rma.mv(yi = yi,
                  V = vi,
                  slab = author,
                  data = meta_df %>% 
                  filter(ocean != "Arctic") %>% 
                  mutate(ocean = fct_rev(ocean)),
                  method = "REML",
                  mods = ~ abs_lat*ocean,
                  random = ~ 1 | paper_num/indiv, 
                  test = "t")
summary(mod_o)
# there's a positive relationship in the pacific but maybe not others?



# Model with just phyto and continuous background_nutrients
mod_bn_p <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = meta_df %>% filter(pp == "Phytoplankton"),
                 method = "REML",
                 mods = ~ lat_center*duration_stand,
                 random = ~ 1 | paper_num/indiv, 
                 test = "t")
summary(mod_bn_p)


## Coral vs Endosymbiont ------------------------------------------------------

# make a df of just coral effect sizes
c_nums <- meta_data %>%
  filter(pp_taxa == "Co")

# just the endosymbiont effect sizes
en_nums <- meta_data %>%
  filter(pp_taxa == "En")

# of the coral papers, which of those papers also recorded an endosymbiont
c_en_nums <- c_nums %>%
  filter(paper_num %in% en_nums$paper_num)

# ok now make a final df of those papers
meta_ce <- meta_data %>%
  filter(paper_num %in% c_en_nums$paper_num) %>%
  filter(pp_taxa == "Co" | pp_taxa == "En") %>%
  droplevels() %>%
  as.data.frame() 

# I can only compare endo vs coral bc there aren't temperate or high N corals
mod_ce <-  rma.mv(yi = yi,
                     V = vi,
                     slab = author,
                     data = meta_ce,
                     method = "REML",
                     mods = ~ pp,
                     random = ~ 1 | paper_num/indiv, 
                     test = "t")

summary(mod_ce)
# write_csv(save_mod_output(mod_ce), "Tables/mod_ce.csv") 


# prep tt model for emmeans
prep_ce_tt <- emmeans_prep(mod_ce_tt) 

# write_csv(save_mod_output(mod_ce_tt), "Tables/mod_ce_tt.csv") 

# What if we just compare coral and endosymbionts
meta_coral <- meta_data %>% 
  filter(pp_taxa == "En" | pp_taxa == "Co" ) %>%
  droplevels()

model_coral <- rma.mv(yi = yi,
                      V = vi,
                      slab = author,
                      data = meta_coral,
                      method = "REML",
                      mods = ~ pp,
                      random = ~ 1 | paper_num/indiv, 
                      test = "t")

summary(model_coral)

prep_mod_co <- emmprep(model_coral) 
pred_mod_co_table <- as.data.frame(emmeans(prep_mod_co, specs="pp",  type="response", tran = NA)) 

ggplot(pred_mod_co_table, aes(emmean, pp, colour = pp)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  theme_white() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)")



## Tropics vs Subtropics vs Temperate ------------------------------------------
# tropics vs subtropics vs temperate

mod_tt2 <- rma.mv(yi = yi,
                  V = vi,
                  slab = author,
                  data = meta_df ,
                  method = "REML",
                  mods = ~ climate*pp + nat_bin2,
                  random = ~ 1 | paper_num/indiv, 
                  test = "t")

summary(mod_tt2)
#plot_summs(mod_tt2)

# write_csv(save_mod_output(mod_tt2), "Tables/mod_tt2.csv") 

# prep model for emmeans
prep_tt2 <- emmeans_prep(mod_tt2) # use function to prep model for emmeans

# pairwise: did pp respond stronger temp vs tropics?
pairs(emmeans(prep_tt2,  ~ climate*pp), simple = "climate") 


### In situ temperature -----------------------------------------------------------
# Model with continuous measured in situ temperature
mod_temp <- rma.mv(yi = yi,
                   V = vi,
                   slab = author,
                   data = meta_df %>% filter(temp_proxy == "measure"),
                   method = "REML",
                   mods = ~ temp_center*pp,
                   random = ~ 1 | paper_num/indiv, 
                   test = "t")

summary(mod_temp)
# write_csv(save_mod_output(mod_temp), "Tables/mod_temp.csv") 

slope_temp <- coef(summary(mod_temp))[2,1] # extract the slope

# Calculate slopes
int_effect_sum_t <- slope_coeff_fun(mod_temp, slope_temp) # FANCY FUNCTION

# Predictions for scatter plot
# Generate values
temp_values <- new_df(temp, meta_df) %>% # Function to create df with dummy vars
  rename(temp = var) %>%
  mutate("temp:ppCoral" = temp*ppCoral,
         "temp:ppEndosymbiont" = temp*ppEndosymbiont,
         "temp:ppPhytoplankton" = temp*ppPhytoplankton,
         "temp:ppSeagrass" = temp*ppSeagrass) %>%
  as.matrix()

# new model with non standardized background nutrients
model_temp_untran <- rma.mv(yi = yi,
                            V = vi,
                            slab = author,
                            data = meta_data,
                            method = "REML",
                            mods = ~ temp + pp + temp:pp,
                            random = ~ 1 | paper_num/indiv, 
                            test = "t")

# use the generated values to plot background nutrients effect for each pp_group, nat_group held to 1
temp_predict <- as.data.frame(predict(object = model_temp_untran,
                                      newmods = temp_values)) %>%
  cbind(temp_values) %>%
  mutate(pp = c(rep("Algae", each = 100), 
                rep("Coral", each = 100), 
                rep("Endosymbiont", each = 100), 
                rep("Phytoplankton", each = 100), 
                rep("Seagrass", each = 100)))

# slope plot
slope_t <- int_effect_plot(int_effect_sum_t) +
  scale_colour_manual(values = rev(pal5)) +
  xlab("Slope") 

# scatter plot
scatter_t <- scatter_fun(temp_predict, meta_df, temp, rev(pal5)) + 
  labs(x = expression("Temperature " ( degree*C))) 

slope_t + scatter_t

### continent? ----
mod_cont <- rma.mv(yi = yi,
                         V = vi,
                         slab = author,
                         data = meta_df %>% mutate(nat_bin2 = fct_rev(nat_bin2)),
                         method = "REML",
                         mods = ~ continent*pp,
                         random = ~ 1 | paper_num/indiv, 
                         test = "t")
summary(mod_cont)

## type of nutrient ----------
mod_poo <- rma.mv(yi = yi,
                  V = vi,
                  slab = author,
                  data = meta_df,
                  method = "REML",
                  mods = ~ prov*tt*pp,
                  random = ~ 1 | paper_num/indiv, 
                  test = "t")

summary(mod_poo)

# prep model for emmeans
prep_poo <- emmeans_prep(mod_poo) 

# pairwise: did pp respond stronger low vs high nutriens?
pairs(emmeans(prep_poo,  ~ prov)) 

pred_poo <- as.data.frame(emmeans(prep_poo, specs="prov", by = c("pp", "tt"), type="response", tran = NA)) 

ggplot(pred_poo, aes(emmean, prov, colour = tt)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  facet_wrap(~pp, ncol= 1, scales = 'free_y') +
  pub_theme() +
  #theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)")


### pee vs poo vs animal -------------------------------------------------------
mod_poo <- rma.mv(yi = yi,
                  V = vi,
                  slab = author,
                  data = meta_df,
                  method = "REML",
                  mods = ~ prov*tt*pp,
                  random = ~ 1 | paper_num/indiv, 
                  test = "t")

summary(mod_poo)

# prep model for emmeans
prep_poo <- emmeans_prep(mod_poo) 

# pairwise: did pp respond stronger low vs high nutrients?
pairs(emmeans(prep_poo,  ~ prov)) 

pred_poo <- as.data.frame(emmeans(prep_poo, specs="prov", by = c("pp", "tt"), type="response", tran = NA)) 

ggplot(pred_poo, aes(emmean, prov, colour = tt)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  facet_wrap(~pp, ncol= 1, scales = 'free_y') +
  pub_theme() +
  #theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)")

# plot 
ugh <- meta_sge %>%
  group_by(pp_taxa, paper_num) %>%
  summarise(mean = mean(yi)) %>%
  full_join(meta_sge, by = c("paper_num", "pp_taxa"))

# plot
ggplot(ugh) +
  geom_point(aes(abs_lat, yi, colour = tt), size = 2, alpha = 0.5) +
  geom_point(aes(abs_lat, mean, colour = pp_taxa), size = 4,  alpha = 0.5)

# how many in each group?
sge_n_bn <- meta_sge %>%
  count(pp_taxa, bn_cat)

sge_n_tt <- meta_sge %>%
  count(pp_taxa, tt)

### Seagrass specific  -----------------------------------

# Do seagrass metrics respond differently?
# metrics: belowground biomass, aboveground biomass, leaf measurement (area, width, length), density of shoots, performance

seagrass <- meta_df %>%
  filter(pp_taxa == "S") %>%
  mutate(type = case_when(measured == "shoot" & response == "length" ~ "leaf meas",
                          measured == "shoot" & response == "biomass" ~ "leaf mass",
                          measured == "shoot" & response == "density" ~ "shoot count",
                          measured == "leaf" & response == "count" ~ "leaf count",
                          measured == "leaf" & response == "length" ~ "leaf meas",
                          measured == "shoot" & response == "photosynthesis" ~ "perf",
                          measured == "aboveground" & response == "biomass" ~ "leaf mass",
                          measured == "leaf" & response == "area" ~ "leaf meas",
                          measured == "belowground" & response == "biomass" ~ "root mass",
                          measured == "leaf" & response == "biomass" ~ "leaf mass",
                          measured == "rhizome" & response == "biomass" ~ "root mass",
                          measured == "root" & response == "biomass" ~ "root mass",
                          measured == "leaf" & response == "width" ~ "leaf meas"))
  
a <- seagrass %>%
  count(measured,response)

mod_seagrass <-  rma.mv(yi = yi,
                       V = vi,
                       slab = author,
                       data = seagrass,
                       method = "REML",
                       mods = ~ -1 +type,
                       random = ~ 1 | paper_num/indiv, 
                       test = "t")

summary(mod_seagrass)
plot_summs(mod_seagrass)

# just compare responses between seagrass species 
mod_seagrass_sp <-  rma.mv(yi = yi,
                        V = vi,
                        slab = author,
                        data = seagrass,
                        method = "REML",
                        mods = ~ -1 +pp_genus,
                        random = ~ 1 | paper_num/indiv, 
                        test = "t")

summary(mod_seagrass_sp)

# what about consumer taxa
mod_seagrass_c <-  rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = seagrass %>% filter(c_group == "Molluscs"),
                           method = "REML",
                           mods = ~ sum_all_N_treat*pp_genus,
                           random = ~ 1 | paper_num/indiv, 
                           test = "t")

summary(mod_seagrass_c)

# what about water column vs sediment nutrient additions
mod_seagrass_sp <-  rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = seagrass,
                           method = "REML",
                           mods = ~ -1 +pp_genus,
                           random = ~ 1 | paper_num/indiv, 
                           test = "t")

summary(mod_seagrass_sp)

seagrass %>% filter(paper_num == "20") %>% ggplot() +
  geom_point(aes(treatment_levels, yi), alpha = 0.5)

# 
pos_neg <- meta_df %>% 
  filter(pp == "Seagrass")%>%
  mutate(sign = ifelse(yi >0, "positive", "negative")) %>% 
  group_by(tt, c_group) %>%
  count(sign) %>%
  ungroup() %>%
  pivot_wider(names_from = sign, values_from = n) %>%
  mutate(percent_pos = 100*positive / (positive + negative))


### Finer primary producer groupings -------------------------------------------
meta_df_a <- meta_df %>%
  filter(pp == "Algae") %>%
  mutate(pp_fine = case_when(pp_fine == "epiphyte_turf" ~ "turf",
                             pp_fine == "epiphyte_calc" ~ "calc",
                             TRUE ~ as.character(pp_fine)
                             ))
mod_pp <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = meta_df_a,
                 method = "REML",
                 mods = ~ -1 + pp_fine,
                 random = ~ 1 | paper_num/indiv, 
                 test = "t")

summary(mod_pp)


# what if we look at lowest possible level?
meta_df_df <- meta_df %>%
  mutate(pp_lowest = coalesce(pp_genus, pp_family),
         pp_lowest = case_when(pp_lowest == "spp." ~ pp_genus,
                               pp_lowest == "sp" ~ pp_genus,
                               TRUE ~ as.character(pp_lowest))) %>%
  drop_na(pp_lowest)

mod_pp2 <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = meta_df,
                 method = "REML",
                 mods = ~ pp_genus,
                 random = ~ 1 | paper_num/indiv, 
                 test = "t")

summary(mod_pp2)

# prep model for emmeans
prep_pp2 <- emmeans_prep(mod_pp2) # use function to prep model for emmeans

# make table
pred_pp_table2 <- as.data.frame(emmeans(prep_pp2, specs="pp_genus", type="response", tran = NA)) %>%
  mutate(y_num = 1:length(df),
         dum = "Overall") %>%
  left_join((meta_df %>%
               group_by(pp_genus) %>%
               summarise(n = n()) %>%
               # this creates the n = label
               mutate(n_lab = sub("^","n = ",n))), by = c("pp_genus")) %>% 
  replace_na(list(n_lab = "n = 0")) %>%
  ungroup() %>% 
  left_join(meta_df_df %>% select(pp_genus, pp, pp_family), by = "pp_genus") %>%
  mutate(pp_family = as.factor(pp_family),
         pp = ifelse(pp == "Endosymbiont", "Endo", as.character(pp)))%>%
  arrange(pp_family, pp_genus) %>%
  mutate(pp_genus = factor(pp_genus, levels = unique(pp_genus)))


# plot
ggplot(pred_pp_table2, aes(x = emmean, 
                           y = pp_genus, 
                           colour = pp_family)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  facet_grid(rows = vars(pp),
             scales = "free_y",
             space = "free_y") +
  pub_theme() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)") +
  # Add text labels inside the plots
  geom_text(aes(x = emmean, label = n), # change to 1 if i add tt
            color = "black", size = 3, vjust = -0.5, hjust = 0) +
  scale_y_discrete(expand = expansion(add = 0.85))
# ggsave("Figures/Fig.sp.png", device = "png", height = 170, width = 168, units = c("mm"), dpi = 600)



# compare cordgrass and mangroves with seagrass
meta_plant <- meta_data %>%
  filter(pp_group == "Plants")

mod_plant <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = meta_plant,
                    method = "REML",
                    mods = ~ -1 + pp,
                    random = ~ 1 | paper_num/indiv, 
                    test = "t")
summary(mod_plant)



### Conservative models---------------------------------------------------------

# just look at conservative pp responses
meta_c <- meta_df %>%
  filter(!(response == "percent cover" & natural == 1)) # removes 26 rows, 19 of which are algae


### effect of fish on different pp ---------------------------------------------

mod_fish <- rma.mv(yi = yi,
                   V = vi,
                   slab = author,
                   data = meta_data %>% filter(c_taxa == "F"),
                   method = "REML",
                   mods = ~ bn_center*pp_group,
                   random = ~ 1 | paper_num/indiv, 
                   test = "t")

summary(mod_fish)

rownames(mod_fish$beta)[1] = "(Intercept)"
prep_fish = qdrg(object = mod_fish) 

pred_fish_table <- as.data.frame(emmeans(prep_fish, specs="bn_center", by = "pp_group", type="response", tran = NA))


# now plot
ggplot(pred_fish_table, aes(emmean, pp_group, colour = pp_group)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  #facet_wrap(~pp_group, ncol= 1) +
  theme_white() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)")



### does metric effect things? -------------------------------------------------
mod_metric <- rma.mv(yi = yi,
                     V = vi,
                     slab = author,
                     data = meta_data,
                     method = "REML",
                     mods = ~ response*bn_cat + natural,
                     random = ~ 1 | paper_num/indiv, 
                     test = "t")

summary(mod_metric)

mod_metric2 <- rma.mv(yi = yi,
                     V = vi,
                     slab = author,
                     data = meta_data,
                     method = "REML",
                     mods = ~ response,
                     random = ~ 1 | paper_num/indiv, 
                     test = "t")

summary(mod_metric2)
# measurements of size ^2 seem to be higher? that's area and percent cover


### seagrass v molluscs? -------------------------------------------------------
model_sg <- rma.mv(yi = yi,
                      V = vi,
                      slab = author,
                      data = meta_data %>% filter(pp_taxa == "S"),
                      method = "REML",
                      mods = ~ c_taxa*tt,
                      random = ~ 1 | paper_num/indiv, 
                      test = "t")

summary(model_sg)


rownames(model_sg$beta)[1] = "(Intercept)"
prep_sg = qdrg(object = model_sg) 

pred_sg_table <- as.data.frame(emmeans(prep_sg, specs="c_taxa", by = "tt", type="response", tran = NA))


# now plot
ggplot(pred_sg_table, aes(emmean, c_taxa, colour = tt)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  #facet_wrap(~measured, ncol= 1, scales = "free") +
  #theme_white() +
  #theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)")



### check for diff in metric between pp groups ---------------------------------
mod_inst_grow <- rma.mv(yi = yi,
                        V = vi,
                        slab = author,
                        data = meta_data,
                        method = "REML",
                        mods = ~ grow_time_end*pp_group,
                        random = ~ 1 | paper_num/indiv, 
                        test = "t")

prep_grow <- emmeans_prep(mod_inst_grow) # use function to prep model for emmeans

pred_mig_table <- as.data.frame(emmeans(prep_grow, specs="grow_time_end", by = "pp_group", type="response", tran = NA))%>%
  left_join(meta_data %>%
              group_by(pp_group, grow_time_end) %>%
              summarise(n = n()) %>%
              ungroup())

# plot
ggplot(pred_mig_table, aes(emmean, grow_time_end, colour = pp_group)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5) +
  facet_wrap(~pp_group, ncol= 1) +
  pub_theme() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = "dashed") +
  ylab("") +
  xlab("Effect size (LnRR)")


# Look at effect size scatter
meta_data %>%
  group_by(natural) %>%
  mutate(mean_effect = mean(vi)) %>%
  ungroup() %>%
  as.data.frame() %>%
  ggplot() +
  geom_jitter(aes(vi, natural, colour = natural)) +
  geom_point(aes(mean_effect, natural)) +
  geom_vline(xintercept=0, color="black", linetype="dashed") +
  theme(legend.position = "null")


### Look for effects of nh4 for CH2 --------------------------------------------
d <- meta_data %>% 
  drop_na(nh4_treatment, nh4_control) %>% 
  filter(yi >0) %>%
  mutate(x_diff = nh4_treatment/nh4_control) %>%
  select(author, title, yi, x_diff, climate_zone, nutrient_levels, nh4_treatment, nh4_control, no3_treatment, no3_control, prov_units) %>%
  filter(x_diff < 1000) %>%
  filter(x_diff >= 1)

### effect size vs experimental design + added N  -----------
mod_nat_added <- rma.mv(yi = yi,
                        V = vi,
                        slab = author,
                        data = meta_added,
                        method = "REML",
                        mods = ~ added_N_log*nat_bin2,
                        random = ~ 1 | paper_num/indiv, 
                        test = "t")

summary(mod_nat_added)

# try to visualize that?
ggplot(meta_added, aes(added_N_log, yi, colour = nat_bin2)) +
  geom_point() +
  geom_smooth(method = lm) +
  geom_vline(xintercept = 0, lty = "dashed") +
  geom_hline(yintercept = 0, lty = "dashed")
# ggsave("Figures/added_N_vs_nat_bin2.png", device = "png", height = 125, width = 173, units = c("mm"), dpi = 600)


# Presentation figures -------------------
# Map for presentations
world2 <- ne_countries(scale = "medium", returnclass = "sf") 

blue2 <- paste("#99B9F1", sep="")


# Map of world
ggplot() +
  geom_sf(data = world2, fill = "white", colour = blue2) +
  # add aesthetic elements
  theme_bw() +
  theme(panel.background = element_rect(fill = blue2),
        panel.grid.major = element_line(color = blue2),
        axis.text = element_text(size = 9, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 9, colour = "black"),
        legend.margin =margin(r=1,l=1,t=1,b=1),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25))


# Map of studies -----
lat_long2 <- meta_data %>%
  distinct(paper_num, lat_dec, long_dec, tt) %>%
  as.data.frame() %>%
  mutate(lat = lat_dec,
         long = long_dec) %>%
  st_as_sf(coords = c("long_dec", "lat_dec")) %>%
  st_set_crs(4326)

# Load simple world map from package "rnaturalearth"
world2 <- ne_countries(scale = "medium", returnclass = "sf")  %>%
  st_set_crs(4326)


# Plot
ggplot() +
  geom_sf(data = world2, fill = "white", color = "black") +
  # add coordinates for each study
  geom_sf(data = lat_long2, pch = 21, 
          colour = "black", fill = pal1, 
          size = 5, alpha = 0.6) +
  theme_black() +
  theme(panel.background = element_rect(fill = blue),
        panel.grid.major = element_line(color = blue),
        axis.text = element_text(size = 9, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 9, colour = "black"),
        legend.margin =margin(r=1,l=1,t=1,b=1),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
  ylim(c(-75, 75)) +
  xlim(c(-162, 164))
# ggsave("Figures/Pub_figs/Fig.2.1.png", width = 16, height = 9, units = "in", dpi = 600)
 

# mess around
d <- meta_df %>%
  filter(ocean == "Pacific" &
           continent == "North America" &
           pp_taxa == "A") %>%
  filter(paper_num != 45)
