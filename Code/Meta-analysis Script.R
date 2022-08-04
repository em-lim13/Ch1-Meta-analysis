# Code to analyze meta data
# Em Lim
# Dec 2021

# Load packages and data -------

library(tidyverse)
library(ggplot2)
library(PNWColors)
library(lmerTest)
library(lme4)
library(visreg)
library(metafor)
library(esc)
library(RColorBrewer)

source(file = "Code/theme_black_code.R")
theme_set(theme_classic())

# Set preferences so I can see all my columns
rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

# mine data from papers using this site
# https://automeris.io/WebPlotDigitizer/


# Load data
# And do a million manipulations......
meta_full <- read_csv("Output/Data/march_meta.csv") %>%
  mutate(sd_treatment = replace_na(sd_treatment, 999), # assigning a dummy to replace
         sd_control = replace_na(sd_control, 999)) %>% # assigning a dummy to replace
  mutate(sd_treatment = ifelse(sd_treatment == 999, se_treatment*sqrt(treatment_n), sd_treatment), # replacing that dummy with the estimated SD
         sd_control = ifelse(sd_control == 999, se_control*sqrt(control_n), sd_control),# replacing that dummy with the estimated SE
         natural = as.factor(natural),
         abs_lat = abs(lat_dec), # creating an absolute value lat
         log_nut = log(background_nut),
         climate2 = factor(ifelse(abs_lat <= 23.45, "Tropics",
                           ifelse(abs_lat <= 35, "Subtropics",
                                  ifelse(abs_lat <= 60, "Temperate",
                                         "Polar"))),
         levels = c("Tropics", "Subtropics", "Temperate", "Polar")), # New column where I don't distinguish between North subtropics/temperates and South subtropics/temperates
         c_activity = ifelse(c_taxa == "Bi" | c_taxa == "F", "Vertebrate",
                        ifelse(c_taxa == "SP" | c_taxa == "B" | c_taxa == "SI", "Sessile Invert", "Mobile Invert")),
         c_activity = factor(c_activity, levels = c("Vertebrate", "Mobile Invert", "Sessile Invert")), # activity levels for consumers
         c_taxa = factor(c_taxa, levels = c("Bi", "F", "Cr", "E", "G", "MI", "B", "SP", "SI")), # reorder c_taxa
         climate_zone = factor(ifelse(lat_dec >= 66.5, "NP",
                                      ifelse(lat_dec >= 35, "NT",
                                      ifelse(lat_dec >= 23.45, "NST",
                                      ifelse(lat_dec >= -23.45, "T",
                                      ifelse(lat_dec >= -35, "SST",
                                      ifelse(lat_dec >= -66.5, "ST",
                                             "SP")))))), 
                                levels = c("ST", "SST", "T", "NST", "NT", "NP"))) # assign climate zones based on a code, makes it easier and less prone to mistakes
# using 23.45 as the cut off bc that's between the cancer and capricorn?
# Another paper set 23.5 - 23.5, and set temperate as 30-40

# Calculate effect sizes -------

##### I think I need to manually exclude the negative means because it's still calculating effect sizes for those
meta_no_neg <- meta_full %>%
  filter(mean_control > 0) %>%
  filter(mean_treatment > 0)

# We lose 5 rows when we remove 0's and 4 when we remove negatives (9 total)
# should I switch to Hedges g????
# Should I add like 0.0000000000001 to all the 0's?
 # NO DO NOT DO THIS IT RUINS EVERYTHING EVERY STATISTICAL ANALYSIS BECOMES P = 0.99

meta <- meta_no_neg %>%
  escalc(measure = "ROM", data = meta_no_neg,
       m1i = mean_treatment, sd1i = sd_treatment, 
       m2i = mean_control, sd2i = sd_control,
       n1i = treatment_n, n2i = control_n)
# ROM for log transformed ratio of means
# log(treatment/control)

# Gather all papers into dataframe where each is counted once -------
# Make a data frame where there's only one row for each paper to graph where the papers are from
meta_reduced <- meta_full %>%
  filter(obs_num == 1) 

# Load papers that I can't use but that mostly pass the criteria
meta_rejected <- read_csv("Output/Data/rejected_papers.csv") %>%
  mutate(abs_lat = abs(lat_dec),
         climate2 = factor(ifelse(abs_lat <= 23.5, "Tropics",
                                  ifelse(abs_lat <= 40, "Subtropics",
                                         ifelse(abs_lat <= 60, "Temperate",
                                                "Polar"))),
                           levels = c("Tropics", "Subtropics", "Temperate", "Polar"))) %>%
  dplyr::select(author, title, year, lat_dec, long_dec, climate_zone, climate2)

meta_all_papers <- meta_reduced %>%
  dplyr::select(author, title, year, lat_dec, long_dec, climate_zone, climate2) %>%
  rbind(meta_rejected) # df with all the papers I'm using and the papers I can't use but that still mostly meet the criteria

# Add chla data to orig meta dataframe  -------
# this file was written in the Chla_extraction.R code file
chla_all_values <- read_csv("Output/Data/all_extracted_chla_updated.csv") 

meta_chl <- meta %>%
  left_join(chla_all_values, by = "indiv")

# What's the max and min chla value?
max <- max(meta_chl$chla_est) # 26.73629
min<- min(meta_chl$chla_est) # 0.05139377

# How big should each bin be?
bin_size <- max/5
upper_bin1 <- bin_size
upper_bin2 <- bin_size*2
upper_bin3 <- bin_size*3
upper_bin4 <- bin_size*4
upper_bin5 <- bin_size*5

# Smaller bins
bin_size_small <- max/3
upper_bin_low <- bin_size_small
upper_bin_mid <- bin_size_small*2
upper_bin_high <- bin_size_small*3


# use this so that it still works even if the max chla value changes

# add chl a estimates to df
# add a binned chl a column
meta_final <- meta_chl %>%
  mutate(chla_bin = ifelse(chla_est <= upper_bin1, "1",
                      ifelse(chla_est <= upper_bin2, "2",
                             ifelse(chla_est <= upper_bin3, "3", 
                                 ifelse(chla_est <= upper_bin4, "4", "5")))),
         chla_log = log(chla_est),
         chla_bin_num = as.numeric(chla_bin),
         chla_bin_small = factor(ifelse(chla_est <= upper_bin_low, "Low",
                                  ifelse(chla_est <= upper_bin_mid, "Medium",
                                         "High")),
                           levels = c("Low", "Medium", "High"))) %>%
  rowwise() %>%
  mutate(sum_all_N_treat = sum(nh4_treatment, no3_treatment, total_nitrogen_treatment, na.rm = TRUE),
         sum_all_N_cont = sum(nh4_control, no3_control, total_nitrogen_control, na.rm = TRUE),
         added_N = na_if((sum_all_N_treat - sum_all_N_cont), 0),
         added_nh4 = na_if((sum(nh4_treatment, na.rm = TRUE) - sum(nh4_control, na.rm = TRUE)), 0)) %>%
  ungroup() %>%
  mutate(pp_group = ifelse(pp_taxa == "CA" | pp_taxa == "FA", "Algae",
                    ifelse(pp_taxa == "En" | pp_taxa == "Co", "Coral",
                    ifelse(pp_taxa == "P", "Phytoplankton",
                    ifelse(pp_taxa == "S", "Plants",
                    "Plants"))))) %>%
  mutate(pp_group = ifelse(pp_group == "Plants" & pp_phyla == "Cnidaria", "Coral",
                           pp_group)) %>%
  mutate(lat_stand = base::scale(abs_lat),
         background_nut_stand = scale(background_nut),
         added_N_stand = scale(added_N),
         temp_stand = scale(temp))
  
# Exploration graphing --------
# Make palette
pal <- pnw_palette("Bay", 6)
pal2 <- pnw_palette("Bay", 4)
pal_star <- pnw_palette("Starfish", 5)

climate_labs <- c("South Temperate", "South Subtropics", "Tropics", "North Subtropics", "North Temperate", "North Polar")

#using meta_reduced so each paper is only counted once
ggplot(meta_all_papers, aes(y = climate_zone, fill = climate_zone)) +
  geom_bar(stat="count") +
  scale_fill_manual(values = rev(pal)) +
  scale_y_discrete(labels = climate_labs) +
  labs( y = "Climate Zone", x = "Number of Papers") + 
  theme_black() + theme(legend.position = "none")

# Using climate2
ggplot(meta_all_papers, aes(y = climate2, fill = climate2)) +
  geom_bar(stat="count") +
  scale_fill_manual(values = rev(pal2)) +
  labs( y = "Climate Zone", x = "Number of Papers") + 
  theme_black() + theme(legend.position = "none")


#ggsave("Output/Figures/climate_zone_bar", device = "png",
#      height = 9, width = 16, dpi = 400)


# GRAPHING CHL A
ggplot(meta_final, aes(y = lat_dec, x = chla_est, colour = climate_zone)) +
  geom_jitter()+
  scale_colour_manual(values = rev(pal), labels = climate_labs) +
  labs(y = "Latitude", x = "Chlorophyll a (mg/m^2)") + 
  theme_black() + 
  theme(legend.position = "none") +
  ylim(c(-73, 73))

#ggsave("Output/Figures/chla_lat", device = "png",
#     height = 9, width = 16, dpi = 400)

ggplot(meta_final, aes(x = abs_lat, y = chla_est, colour = climate_zone)) +
  geom_jitter()+
  scale_colour_manual(values = rev(pal), labels = climate_labs) +
  labs(x = "Latitude", y = "Chlorophyll a (mg/m^2)") + 
  theme_black() + 
  theme(legend.position = "none") 


# take a look at the distribution of primary producers across latitudes

# this df filters out any rows that can't be used for added_N
meta_smol <- meta_final %>%
  filter(prov_units == "umol/L") %>%
  filter(added_N < 10) %>% # Filter out the super high ammonium estimates that are dragging the line down?
  drop_na(background_nut_stand) %>%
  drop_na(temp_stand) %>%
  drop_na(added_N_stand) 


ggplot(meta_smol, aes(pp_taxa, climate2, colour = climate2)) +
  geom_jitter(size = 4) +
  theme_black() +
  theme(legend.position = "none") +
  labs(x = "Primary Producer", y = "Climate Zone") +
  scale_x_discrete(labels = c("Calcareous algae", "Coral", "Endosymbiont", "Fleshy Algae", "Mangrove", "Phytoplankton", "Seagrass")) +
  theme(axis.text.x = element_text(angle = 90))

# Algae not in tropics, and Coral not in Temperates. Phytoplankton and Plants are spreak
  
# take a look at the distribution of consumers across latitudes
ggplot(meta_final, aes(c_taxa, climate2, colour = c_taxa)) +
  geom_jitter(size = 4) +
  theme_black() +
  theme(legend.position = "none") +
  labs(x = "Consumer", y = "Climate Zone") +
  scale_x_discrete(labels = c("Birds", "Fish", "Crustaceans", "Echinoderms", "Gastropods", "Misc mobile inverts", "Bivalves", "Sponges", "Misc sessile inverts")) +
  theme(axis.text.x = element_text(angle = 30))

# look at consumer groups
ggplot(meta_final, aes(c_activity, climate2, colour = c_activity)) +
  geom_jitter(size = 4) +
  theme_black() +
  theme(legend.position = "none") +
  labs(x = "Consumer Activity Level", y = "Climate Zone")

# Look at consumers vs pp??
ggplot(meta_final, aes(c_taxa, pp_taxa, colour = climate2)) +
  geom_jitter(size = 4) +
  theme_bw() +
  labs(x = "Consumer Taxa", y = "Primary Producer Taxa", colour = "Climate") +
  scale_x_discrete(labels = c("Birds", "Fish", "Crustaceans", "Echinoderms", "Gastropods", "Misc mobile inverts", "Bivalves", "Sponges", "Misc sessile inverts"))+
  theme(axis.text.x = element_text(angle = 30)) +
  scale_y_discrete(labels = c("Calcareous algae", "Coral", "Endosymbiont", "Fleshy Algae", "Mangrove", "Phytoplankton", "Seagrass")) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)) +
  scale_colour_manual(values = rev(pal2))
  
# a lot of bivalve papers coming out of temperate regions

# Looking at background nutrients and temperature and provisioning

ggplot(meta_final, aes(abs_lat, background_nut)) +
  geom_jitter(aes(colour = climate2)) +
  geom_smooth(method = lm) +
  scale_colour_manual(values = rev(pal2)) +
  theme_black() +
  labs(x = "Absolute Latitude", y = "Background nutrients (umol/L)") +
  theme(legend.position = "none")

#ggsave("Output/Figures/abs_lat_nutrients", device = "png",
#     height = 9, width = 16, dpi = 400)

ggplot(meta_final, aes(abs_lat, temp)) +
  geom_point() +
  geom_smooth(method = lm)

ggplot(meta_final, aes(background_nut, chla_est)) +
  geom_point(aes(colour = climate2)) +
  geom_smooth(method = lm)

# Multilevel models ------
# Link to book with help:
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/metareg.html

# There is a special way through which we can tell rma.mv to assume such nested random effects. We do this using a slash (/) to separate the higher- and lower-level grouping variable. To the left of /, we put in the level 3 (cluster) variable. To the right, we insert the lower-order variable nested in the larger cluster. Therefore, the general structure of the formula looks like this: ~ 1 | cluster/effects_within_cluster

# Look at Test of Moderators to determine whether there are significant differences between the subgroups
# The Model Results are printed within a meta-regression framework. This means that we cannot directly extract the estimates in order to obtain the pooled effect sizes within subgroups.


# One model to rule them all ------
# If goal is to explain patterns, and sample size is big enough to include all the predictors, then just throw them all in the same model. 
# Standardize the numeric predictors, then compare standardized effect sizes
# AIC is good for when you can’t put everything in the same model, I might not be able to put temperature and lat in the same model
# But I lose like 100 papers when I drop all the NA’s for provisioned nutrients
# And even when I drop added_N, which increases the sample size, then the # main effect is basically the same
# Added_N doesn't seem to matter as a mod? No matter how I run it, it has a teeeeeny effect and it's not significant

model_all <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_final,
                          method = "REML",
                          mods = ~ lat_stand + background_nut_stand + added_N_stand + temp_stand + natural + pp_group + lat_stand:background_nut_stand + added_N_stand:background_nut_stand +  temp_stand:background_nut_stand + temp_stand:added_N_stand,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_all)
vif(model_all) # Lat and temp are correlated
# intercept tells me what the value is at the natural 1 pp group algae and mean of all other variables
# looking at the full data (took out all added_N terms) the coefficients are basically hte same
# so we can look at the reduced data and it's saying something very similar to the full data

# In general, all three ways I’ve tried to do these stats give a similar answer. Identity of the pp matters (algae has a stronger response than others), and background nutrients is weak/moderate positive effect. This makes both # Hannah and I feel better. All methods point to context dependence!

# Correlation between variables
cor(cbind(meta_final$lat_stand, meta_final$background_nut_stand, meta_final$added_N_stand, meta_final$temp_stand), use = "complete.obs")
# so temp, lat, and background_nutrients are correlated (|value| < 0.5)
# So maybe we should go back to AIC because we shouldn't have these three in the same model
# BC they're correlated, there might be a generally negative relationship between lat/temp/back_nut and effect but within that slope the model things back_nut has a positive effect

#### Three separate models ------
# 1) everything in full model - temp and lat
# 2) everything in full model - temp and back_nut
# 3) everything in full model - back_nut and lat

# Then compete them using AIC!

# rational for AIC is that we can't put all three in the same model bc they mess up each other
# AIC will tell us which of the three is best

#Subset each plant type to rerun the model_all with each plant type? Or just to do algae by itself and then the rest by themselves. 
# -might have to run each mod by itself bc low sample size
# 

# Except I'm going to compete the "full" model with lat vs "full" with temp vs "full" with background_nutrients

# Want all models to use the same data!
# this df filters out any rows that can't be used for added_N
meta_smol <- meta_final %>%
  filter(prov_units == "umol/L") %>%
  filter(added_N < 10) %>% # Filter out the super high ammonium estimates that are dragging the line down?
  drop_na(background_nut_stand) %>%
  drop_na(temp_stand) %>%
  drop_na(added_N_stand) 

# this df only makes sure to standardize the df between the three correlated variables, and is for use in models without added_N
meta_smolish <- meta_final %>%
  drop_na(background_nut_stand) %>%
  drop_na(temp_stand) 

model_full_lat <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = meta_smolish,
                    method = "REML",
                    mods = ~ lat_stand  + natural + pp_group,
                    random = ~ 1 | paper_num/indiv, 
                    test = "t")

summary(model_full_lat)

# Full model with background_nutrients
model_full_back <- rma.mv(yi = yi,
                                V = vi,
                                slab = author,
                                data = meta_smolish,
                                method = "REML",
                                mods = ~ background_nut_stand  + natural + pp_group,
                                random = ~ 1 | paper_num/indiv, 
                                test = "t")

summary(model_full_back)


# AIC(model_full_back, model_full_back_added)  # make sure to run this with meta_smol for both models
# They're exactly the same. adding added_N doesn't improve the model at all so I want to drop it to increase sample size..?
# The model with added_N and a smaller sample size comes out looking EXACTLY like the model without added_N and a larger sample size. Let's keep added_N for now
# The only thing that changes is that coefficient for back_nut becomes signif when we get rid of added_N and increase sample size


# Full model with temperature
model_full_temp <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_smolish,
                          method = "REML",
                          mods = ~ temp_stand + natural + pp_group,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_full_temp)

# Null model
model_null <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_smolish,
                          method = "REML",
                          mods = ~ 1,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_null)

# Model with just pp and nat
model_cat <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_smolish,
                          method = "REML",
                          mods = ~ pp_group + natural,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_cat)

# 
AIC_tab <- AIC(model_full_lat, model_full_back, model_full_temp, model_null, model_cat) %>%
  rownames_to_column() %>%
  mutate(best = min(AIC),
         delta = AIC - best)

# So the model with background nutrients is the best model according to AIC but there's still moderate support for all three of lat, temp, and background nutrients
# Still the best model even if I get rid of added_N which increases the sample size
# Then present the coefficient plot of the best model, but include the other two in the supp
# And we still get an estimate for the effect size for everything
# Since none of the models are better than the other, I can just chose the one I'm most interested in and throw the other two in the supp
# But.... I'm interested in both lat and background nut
# Tried running interaction models between added_N and lat/temp/back_nut and the interaction didn't change anything so I'll leave it out for simplicity. 

# I'm tempted to keep the background
# Sort of tempted to drop added_N because it restricts the sample size so much and it literally doesn't change the AIC or effects at all and has a teeeeny effect size itself...

# AIC tells you which collection of variables best explain the data
# And then looking at the p-values of those variables tells you whether ___
# Hannah advocates for just looking at the effect size (coefficient)

# Three models with added_N -------
# Lat
model_full_lat_added <- rma.mv(yi = yi,
                         V = vi,
                         slab = author,
                         data = meta_smol,
                         method = "REML",
                         mods = ~ lat_stand  + added_N + natural + pp_group,
                         random = ~ 1 | paper_num/indiv, 
                         test = "t")

summary(model_full_lat_added)

# Full model with background_nutrients
model_full_back_added <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_smol,
                          method = "REML",
                          mods = ~ background_nut_stand + added_N + natural + pp_group,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_full_back_added)


# Full model with temperature
model_full_temp_added <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_smol,
                          method = "REML",
                          mods = ~ temp_stand + added_N + natural + pp_group,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_full_temp_added)

# 
AIC_tab <- AIC(model_full_lat, model_full_back, model_full_temp, model_null, model_cat) %>%
  rownames_to_column() %>%
  mutate(best = min(AIC),
         delta = AIC - best)

# 3 plots: --------
# Plot with pp_groups coefficient + CI
# Plot with natural coefficient + CI
# Plot with coefficient + CI of each continuous variable

# Maybe some help here?
# https://www.metafor-project.org/doku.php/tips:meta_regression_with_log_rr

# this is the final model i swear to god please let this be the final model
model_full_back <- rma.mv(yi = yi,
                                V = vi,
                                slab = author,
                                data = meta_smolish,
                                method = "REML",
                                mods = ~ background_nut_stand + natural + pp_group,
                                random = ~ 1 | paper_num/indiv, 
                                test = "t")

summary(model_full_back)

# pred <- as.data.frame(predict(model_full_back_added, newmods=c(1,0), transf=exp, digits=2))
# coef <- coef(summary(model_full_back_added))

# Plot PP first -------
# Graph primary producer
# How many papers in each bin?
pp_n <- meta_smolish %>%
  group_by(pp_group) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:4)

papers <- meta_smolish %>%
  group_by(paper_num) %>%
  summarise()

# Pull coefficient estimates so I can graph them
intercept_pp <- coef(summary(model_full_back))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes

pp_effect_sum <- coef(summary(model_full_back)) %>%
  slice(-(2:6)) %>%
  mutate(pp_group = factor(c("Algae", "Coral", "Phytoplankton", "Plants")), 
         effect = ifelse(pp_group != "Algae", estimate + intercept_pp,
                         intercept_pp), # add intercept to the estimated coefficient
         ci.lb_adjust = effect - (1.96*se), #adjust the CI based on the new estimate
         ci.up_adjust = effect + (1.96*se)) %>% 
  left_join(pp_n, by = "pp_group")


# GRAPH THIS?
pp_labs <- c("Algae", "Coral", "Phytoplankton", "Plants")

pp_n_labs <- c("n = 25", "n = 47", "n = 53", "n = 47")

ggplot(pp_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb_adjust, xmax = ci.up_adjust, colour = pp_group)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept = 0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size", y = " ") +
  scale_y_continuous(breaks = 1:length(pp_labs),
                     labels = pp_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(pp_n_labs),
                                         labels = pp_n_labs)) +
  scale_colour_manual(values = rev(pal_star)) +
  theme_black()+ 
  theme(legend.position = "none")

# ggsave("Output/Figures/pres_pp.png", device = "png",
#      height = 9, width = 16, dpi = 400)


# Graph experimental design -----
# How many papers in each bin?
nat_n <- meta_smolish %>%
  group_by(natural) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:5,
         nat_bin = factor(c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm"),
                          levels = c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")))

# Pull coefficient estimates so I can graph them
intercept_nat <- coef(summary(model_full_back))[1,1] # extract the intercept

# When you add intercept to the other estimates, the p-value asks whether the CI’s intercept the intercept, not whether they intercept 0

# Can’t just add intercept to CI, works for the estimates but not the CI
# Try using emmeans for this

# edit bin labels
nat_effect_sum <- coef(summary(model_full_back)) %>%
  slice(-(2)) %>%
  slice(-(6:8)) %>%
  mutate(nat_bin = factor(c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")), # add intercept to the estimated coefficient
         effect = ifelse(nat_bin != "Field: unstocked", estimate + intercept_nat,
                         intercept_nat), # add intercept to the estimated coefficient
         ci.lb_adjust = effect - (1.96*se), #adjust the CI based on the new estimate
         ci.up_adjust = effect + (1.96*se)) %>% 
  left_join(nat_n, by = "nat_bin")

# GRAPH THIS?
nat_labs <- c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")

nat_n_labs <- c("n = 74", "n = 14", "n = 21", "n = 9", "n = 54")

ggplot(nat_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb_adjust, xmax = ci.up_adjust, colour = nat_bin)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept = 0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size", y = " ") +
  scale_y_continuous(breaks = 1:length(nat_labs),
                     labels = nat_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(nat_n_labs),
                                         labels = nat_n_labs)) +
  scale_colour_manual(values = rev(pal_star)) +
  theme_black()+ 
  theme(legend.position = "none")

# ggsave("Output/Figures/pres_nat.png", device = "png",
#      height = 9, width = 16, dpi = 400)


# Graph continuous variables -----

# Pull coefficient estimates so I can graph them
lat_effect_sum <- coef(summary(model_full_lat)) %>%
  slice(-(3:9)) %>%
  slice(-1) %>%
  rownames_to_column() 

temp_effect_sum <- coef(summary(model_full_temp)) %>%
  slice(-(3:9)) %>%
  slice(-1) %>%
  rownames_to_column() %>%
  rbind(lat_effect_sum)

cont_effect_sum <- coef(summary(model_full_back)) %>%
  slice(-(3:9)) %>%
  slice(-1)%>%
  rownames_to_column() %>%
  rbind(temp_effect_sum) %>%
  mutate(dummy = 1)


# GRAPH THIS?
cont_labs <- c("Background nutrients", "Absolute latitude", "Temperature")
names(cont_labs) <- c("background_nut_stand", "lat_stand", "temp_stand")

ggplot(cont_effect_sum, aes(y = dummy, x = estimate, xmin = ci.lb, xmax = ci.ub, colour = rowname)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Standardized effect size", y = " ")  +
  theme_black()+ 
  theme(legend.position = "none") +
  scale_y_discrete(labels = cont_labs) +
  scale_colour_manual(values = rev(pal_star)) +
  facet_wrap(~ rowname, nrow = 3, labeller = labeller(rowname = cont_labs))

# ggsave("Output/Figures/pres_cont.png", device = "png",
#      height = 9, width = 16, dpi = 400)

# Independence and Bias -------

funnel.rma(model_lat)
funnel.rma(model_back_nutrients)
funnel.rma(model_added_N)
funnel.rma(model_temp)
funnel.rma(model_exp_design)
funnel.rma(model_act)
funnel.rma(model_pp)

# full model 
funnel.rma(model_full_back)

# A bit of a right skew

fsn(yi = yi, vi = vi, data = meta_smol, type="Rosenthal") #file drawer analysis (missing studies have an effect size = 0)
# Fail-safe N: 13,962
# That seems like quite a large number?????

# Check for influential studies? ------

# make a df where each row is a paper, and num_obs is the number of observations that paper contributed
max_obs <- meta_smolish %>%
  group_by(paper_num) %>%
  dplyr::summarise(num_obs = as.numeric(n()))

hist(max_obs$num_obs)

total_num_obs <- sum(max_obs$num_obs)

# there's 172 data points
# so 10% of that would be 17.2 observations per paper
# so that would mean that paper # 17 and #23 and #29 have too many

# Ok so let's remove those papers one at a time and rerun the stats!

# Removing #29 makes the effect of background_nutrients slightly smaller. 

# # 17 and #23 and #29 are the three I need to check
meta_smolish_drop <- meta_smolish %>%
  filter(paper_num != 29) 

# Rerun with lat model
model_full_lat_drop <- rma.mv(yi = yi,
                               V = vi,
                               slab = author,
                               data = meta_smolish_drop,
                               method = "REML",
                               mods = ~ lat_stand  + natural + pp_group,
                               random = ~ 1 | paper_num/indiv, 
                               test = "t")

summary(model_full_lat_drop)

# Original model output: 
# Test of Moderators (coefficients 2:9):
#   F(df1 = 8, df2 = 163) = 1.7103, p-val = 0.0995

# Model Results:
#                         estimate  se       tval   df    pval    ci.lb    ci.ub 
# intrcpt                  0.5936  0.1608   3.6922  163  0.0003   0.2761   0.9110  *** 
# lat_stand                0.0151  0.0784   0.1919  163  0.8480  -0.1398   0.1699      
# natural2                -0.2354  0.2317  -1.0161  163  0.3111  -0.6930   0.2221      
# natural3                 0.3407  0.2191   1.5548  163  0.1219  -0.0920   0.7735      
# natural4                 0.5877  0.3145   1.8688  163  0.0634  -0.0333   1.2087    . 
# natural5                 0.1651  0.2037   0.8107  163  0.4187  -0.2371   0.5673      
# pp_groupCoral           -0.2233  0.2134  -1.0463  163  0.2970  -0.6448   0.1981      
# pp_groupPhytoplankton   -0.5246  0.2443  -2.1472  163  0.0333  -1.0070  -0.0422    * 
# pp_groupPlants          -0.5219  0.1999  -2.6113  163  0.0099  -0.9166  -0.1273   ** 

# Drop #17 effect of lat goes from 0.0151 to -0.0248. p-val still super high
            # basically there's no effect of lat, 0.01 and -0.02 are reallyyyy straddling 0
# Drop #23 nothing happens
# Drop #29 effect of lat goes from 0.0151 to -0.0218, p-value still super high
            # again, it's still straddling zero pretty hard. it's still basically no effect

# Rerun with background nut model
model_full_back_drop <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_smolish_drop,
                          method = "REML",
                          mods = ~ background_nut_stand  + natural + pp_group,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_full_back_drop)

# Without losing anything: 
# Test of Moderators (coefficients 2:9):
#   F(df1 = 8, df2 = 163) = 2.1250, p-val = 0.0362

# Model Results:
#                         estimate     se     tval   df    pval    ci.lb    ci.ub 
# intrcpt                  0.5881  0.1569   3.7477  163  0.0002   0.2782   0.8980  *** 
# background_nut_stand     0.1274  0.0739   1.7238  163  0.0866  -0.0185   0.2733    . 
# natural2                -0.2114  0.2265  -0.9333  163  0.3521  -0.6586   0.2359      
# natural3                 0.4273  0.2128   2.0078  163  0.0463   0.0071   0.8476    * 
# natural4                 0.6215  0.3111   1.9976  163  0.0474   0.0072   1.2358    * 
# natural5                 0.0940  0.2029   0.4631  163  0.6439  -0.3067   0.4947      
# pp_groupCoral           -0.1564  0.2009  -0.7788  163  0.4373  -0.5531   0.2402      
# pp_groupPhytoplankton   -0.5010  0.2367  -2.1168  163  0.0358  -0.9683  -0.0336    * 
# pp_groupPlants          -0.5510  0.1949  -2.8274  163  0.0053  -0.9358  -0.1662   ** 

# Dropping #17 natural2 coefficient goes from -0.2114 to 0.3155, nothing else changes
# Dropping #23 natural5 coefficient goes from 0.0940 to -0.0079, nothing else changes
              # not a significant difference between nat5 and nat1 either way
# Dropping #29 changes back_nut coefficient from 0.1274 to 0.0670
              # and the p-value goes from 0.0866 to 0.3947
              # everything else is the same

# Full model with temperature
model_full_temp_drop <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_smolish_drop,
                          method = "REML",
                          mods = ~ temp_stand + natural + pp_group,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_full_temp_drop)

# Without losing anything: 
#Test of Moderators (coefficients 2:9):
#  F(df1 = 8, df2 = 163) = 1.7926, p-val = 0.0819

#Model Results:
#                        estimate    se     tval   df    pval    ci.lb    ci.ub 
#intrcpt                  0.5804  0.1606   3.6134  163  0.0004   0.2632   0.8976  *** 
#temp_stand              -0.0554  0.0721  -0.7681  163  0.4436  -0.1979   0.0870      
#natural2                -0.2180  0.2294  -0.9501  163  0.3435  -0.6711   0.2351      
#natural3                 0.3643  0.2144   1.6992  163  0.0912  -0.0590   0.7877    . 
#natural4                 0.6081  0.3142   1.9354  163  0.0547  -0.0123   1.2285    . 
#natural5                 0.1591  0.2021   0.7873  163  0.4323  -0.2399   0.5581      
#pp_groupCoral           -0.1873  0.2089  -0.8968  163  0.3711  -0.5998   0.2251      
#pp_groupPhytoplankton   -0.5221  0.2402  -2.1738  163  0.0312  -0.9963  -0.0478    * 
#pp_groupPlants          -0.5218  0.1967  -2.6520  163  0.0088  -0.9103  -0.1333   ** 

# Dropping #17 natural2 coefficient goes from -0.2180 to 0.3455
     # Coral coefficient goes from -0.1873 to -0.4487 and p-val goes from 0.3711 to 0.0596
# Dropping #23 doesn't change anything
# Dropping #29 doesn't change anything

# " " -----------

# CODE GRAVEYARD ---------

# Chlorophyll a multilevel models -----
# Originally we were using chl a as a proxy for background nutrients. But now I've actually found estimates for background nutrients, and determined that background nutrients and chl a aren't actually that correlated. Which leaves me wondering what chl a is actually a predictor of, knee-jerk would be productivity but we know tropical regions are super productive and have very low chl a
# So maybe we don't look at chl a anymore........

# Ok there are multiple ways to look at chl a data
# We could look at the untransformed chl a, log transformed chl a, numeric binned chl a, or binned categorical chl a
# The only time we see a siginificant effect is when we look at categorical, and that model is probably overfitting the bins with very few samples

# Untransformed chl a estimates as a linear model
model_chla_est <- rma.mv(yi = yi,
                         V = vi,
                         slab = author,
                         data = meta_final,
                         method = "REML",
                         mods = ~ chla_est,
                         random = ~ 1 | paper_num/indiv,
                         test = "t")

summary(model_chla_est)
# p value is 0.5450

# Transformed log chl a data
# Do this when you know that a predictor wouldn't increase the response forever (non linear)
model_chla_log <- rma.mv(yi = yi,
                         V = vi,
                         slab = author,
                         data = meta_final,
                         method = "REML",
                         mods = ~ chla_log,
                         random = ~ 1 | paper_num/indiv,
                         test = "t")
summary(model_chla_log)
# p value is 0.9213

# Binned chl a as numeric
model_chla_paper_numeric <- rma.mv(yi = yi,
                                   V = vi,
                                   slab = author,
                                   data = meta_final,
                                   method = "REML",
                                   mods = ~ chla_bin_num,
                                   random = ~ 1 | paper_num/indiv,
                                   test = "t")
summary(model_chla_paper_numeric)
# p value is 0.5423

# Binned categorical chl a
model_chla_paper_bin <- rma.mv(yi = yi,
                               V = vi,
                               slab = author,
                               data = meta_final,
                               method = "REML",
                               mods = ~ chla_bin,
                               random = ~ 1 | paper_num/indiv, 
                               test = "t")

summary(model_chla_paper_bin)
# p value is 0.0111


AIC(model_chla_paper_numeric, model_chla_paper_bin, model_chla_log, model_chla_est)
# ok so interestingly the bins are the lowest AIC
# sample sizes are super uneven
# This is probably why the AIC is lower in the bins
# Probably overfitting the bins with very few samples
# Hannah says that as long as I justify what I did i might be ok
# But she would be sketched out if she read that i binned
# Hannah suggest going with the continuous chla est
# That's her gut instinct


# Basically Hannah says that the way it's currently coded, the model is just comparing each bin, but doesn't acknowledge that bins 1 and 2 are more similar than 1 and 5
# Could go the ordinal route, might be conceptually ok but logistically challenging to interpret
# links for that:
# https://www.researchgate.net/publication/318745880_Regression_with_Ordered_Predictors_via_Ordinal_Smoothing_Splines
# https://cran.r-project.org/web/packages/ordPens/ordPens.pdf

# OK I want to use chla_est

# Graph chl a est------
chla_predict <- as.data.frame(predict(model_chla_est,
                                      addx=TRUE)) %>%
  tibble::rownames_to_column("author") %>%
  rename(chla_est = X.chla_est) 

# set seed because otherwise I'll get different colours each time
set.seed(1235)

# plot predicted values as line, jitter observed effect sizes
ggplot() +
  geom_line(data = chla_predict, aes(chla_est, pred), colour = "white", size = 1.5) +
  geom_jitter(data = meta_final, aes(chla_est, yi, colour = author)) +
  geom_ribbon(data = chla_predict, aes(chla_est, pred, ymin = ci.lb, ymax = ci.ub), fill = "white", alpha = .1) +
  theme_classic() +
  labs(x = expression(paste(Chlorophyll~a~(mg/m^{2}))), y = "Effect Size") +
  scale_color_manual(values = sample(col_vector, n)) +
  geom_hline(yintercept=0, color="white", linetype="dashed") +
  theme_black() +
  theme(legend.position="none")


#ggsave("Output/Figures/chla_est_effect", device = "png",
#     height = 9, width = 16, dpi = 400)


# Chl a * latitude interaction multi level model -------
# this might be sketchy be chl a is sort of correlated with lat??????
model_multi <- rma.mv(yi = yi,
                      V = vi,
                      slab = author,
                      data = meta_final,
                      method = "REML",
                      mods = ~ chla_est * abs_lat,
                      random = ~ 1 | paper_num/indiv,
                      test = "t")

summary(model_multi)
# Ok test of mods is p value = 0.0194
# chla_est is p value = 0.0039
# lat is p value = 0.8757
# and interaction is p value = 0.0027

# check for colinearity
vif(model_multi) # run this on the model without an interaction
# value of 1 means they're not correlated with any other variables
# higher value = greater correlation
# values of more than 4 or 5 are regarded as moderate - high
# values of 10+ are very high

# OKKKKKKK so basically there's an interaction between lat and chl a ????
# what does that mean
# negative interaction tells us that the two variables have an antagonistic effect – increasing both variables simultaneously results in a smaller effect size than if we added the main effects of each together

# so the further you get from the poles, and the higher the chl a is, the less effect pee has?

# LOG CHLA?
model_multi_log <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_final,
                          method = "REML",
                          mods = ~ chla_log * abs_lat,
                          random = ~ 1 | paper_num/indiv,
                          test = "t")

summary(model_multi_log)

AIC(model_multi, model_multi_log)
#OOOOOKAY so the untransformed model is slightly better

# Graph interaction model -------

# make own df of lat values
mean_lat_chla <- meta_final %>%
  group_by(climate2) %>%
  summarise(mean = mean(abs_lat),
            min_chla = min(chla_est),
            max_chla = max(chla_est))

# what are the min and max chl a values for each zome


lat_values <- data.frame(chla_est = c(seq(from = 0.09748242, to = 1.522176, length.out = 20),
                                      seq(from = 0.05139377, to = 6.724741, length.out = 20),
                                      seq(from = 0.85469812, to = 26.736286, length.out = 20),
                                      seq(from = 0.44598898, to = 0.445989, length.out = 20)),
                         abs_lat = c(rep(17.42045, times = 20), 
                                     rep(26.54635, times = 20), 
                                     rep(48.00235, times = 20), 
                                     rep(74.21667, times = 20))) %>%
  mutate(`chla_est:abs_lat` = chla_est*abs_lat) %>%
  as.matrix()

# change the sequence for chl a to represent the chl a values we actually see at these latitudes (similar notation to how I did latitude with c() instead of rep())
# and change lat values for each bin to better reflect the avg lat in that bin

# use predict to tell it which values you want predictions for
lat_chla_predict <- as.data.frame(predict(object = model_multi,
                                          newmods = lat_values)) %>%
  cbind(lat_values) %>%
  mutate(climate = factor(c(rep("Tropics", times = 20), 
                            rep("Subtropics", times = 20), 
                            rep("Temperate", times = 20), 
                            rep("Polar", times = 20)),
                          levels = c("Tropics", "Subtropics", "Temperate", "Polar")))

# plot predicted values as line, jitter observed effect sizes
ggplot() +
  geom_line(data = lat_chla_predict,  aes(chla_est, pred, colour = climate), size = 1.7) +
  geom_point(data = meta_final, aes(chla_est, yi, colour = climate2), size = 2.5) +
  geom_ribbon(data = lat_chla_predict, aes(chla_est, pred, ymin = ci.lb, ymax = ci.ub, fill = climate), alpha = .2)+
  labs(x = expression(paste(Chlorophyll~a~(mg/m^{2}))), y = "Effect Size", colour = "Climate", fill = "Climate") +
  geom_hline(yintercept=0, color="white", linetype="dashed") +
  theme_black() +
  xlim(c(0, 21.5)) + ylim(c(-1, 3.2)) +
  scale_colour_manual(values = rev(pal2), labels = c("Tropics n = 37", 
                                                     "Subtropics n = 64",
                                                     "Temperate n = 70", 
                                                     "Polar n = 2")) +
  scale_fill_manual(values = rev(pal2), labels = c("Tropics n = 37", 
                                                   "Subtropics n = 64",
                                                   "Temperate n = 70", 
                                                   "Polar n = 2")) 

#ggsave("Output/Figures/interaction_effects", device = "png",
#     height = 9, width = 16, dpi = 400)

# Just plot tropics
trop_predict <- lat_chla_predict %>%
  filter(climate == "Tropics")

trop_meta <- meta_final %>%
  filter(climate2 == "Tropics")

# plot just tropics
ggplot() +
  geom_line(data = trop_predict,  aes(chla_est, pred, colour = climate), size = 1.7) +
  geom_point(data = trop_meta, aes(chla_est, yi, colour = climate2), size = 2.5) +
  geom_ribbon(data = trop_predict, aes(chla_est, pred, ymin = ci.lb, ymax = ci.ub, fill = climate), alpha = .2) +
  labs(x = expression(paste(Chlorophyll~a~(mg/m^{2}))), y = "Effect Size", colour = "Climate", fill = "Climate") +
  geom_hline(yintercept=0, color="white", linetype="dashed") +
  theme_black() +
  xlim(c(0, 21.5)) + ylim(c(-1, 3.2))+
  scale_colour_manual(values = rev(pal2), labels = c("Tropics n = 37", 
                                                     "Subtropics n = 64",
                                                     "Temperate n = 70", 
                                                     "Polar n = 2")) +
  scale_fill_manual(values = rev(pal2), labels = c("Tropics n = 37", 
                                                   "Subtropics n = 64",
                                                   "Temperate n = 70", 
                                                   "Polar n = 2")) 

#ggsave("Output/Figures/interaction_effects_trop", device = "png",
#     height = 9, width = 16, dpi = 400)

# plot just tropics and subtropics
trop_sub_predict <- lat_chla_predict %>%
  filter(climate %in% c("Tropics", "Subtropics"))

trop_sub_meta <- meta_final %>%
  filter(climate2 %in% c("Tropics", "Subtropics"))

# which studies have high chl a in subtropics?
sub_chla <- trop_sub_meta %>%
  filter(chla_est >5)

# plot just tropics and subtropics
ggplot() +
  geom_line(data = trop_sub_predict,  aes(chla_est, pred, colour = climate), size = 1.7) +
  geom_point(data = trop_sub_meta, aes(chla_est, yi, colour = climate2), size = 2.5) +
  geom_ribbon(data = trop_sub_predict, aes(chla_est, pred, ymin = ci.lb, ymax = ci.ub, fill = climate), alpha = .2) +
  labs(x = expression(paste(Chlorophyll~a~(mg/m^{2}))), y = "Effect Size", colour = "Climate", fill = "Climate") +
  geom_hline(yintercept=0, color="white", linetype="dashed") +
  theme_black() +
  xlim(c(0, 21.5)) + ylim(c(-1, 3.2)) +
  scale_colour_manual(values = rev(pal2), labels = c("Tropics n = 37", 
                                                     "Subtropics n = 64",
                                                     "Temperate n = 70", 
                                                     "Polar n = 2")) +
  scale_fill_manual(values = rev(pal2), labels = c("Tropics n = 37", 
                                                   "Subtropics n = 64",
                                                   "Temperate n = 70", 
                                                   "Polar n = 2")) 

#ggsave("Output/Figures/interaction_effects_trop_sub", device = "png",
#     height = 9, width = 16, dpi = 400)

# plot tropics, subtropics, and temperate
trop_sub_temp_predict <- lat_chla_predict %>%
  filter(climate %in% c("Tropics", "Subtropics", "Temperate"))

trop_sub_temp_meta <- meta_final %>%
  filter(climate2 %in% c("Tropics", "Subtropics", "Temperate"))

# plot just tropics and subtropics
ggplot() +
  geom_line(data = trop_sub_temp_predict,  aes(chla_est, pred, colour = climate), size = 1.7) +
  geom_point(data = trop_sub_temp_meta, aes(chla_est, yi, colour = climate2), size = 2.5) +
  geom_ribbon(data = trop_sub_temp_predict, aes(chla_est, pred, ymin = ci.lb, ymax = ci.ub, fill = climate), alpha = .2) +
  labs(x = expression(paste(Chlorophyll~a~(mg/m^{2}))), y = "Effect Size", colour = "Climate", fill = "Climate") +
  geom_hline(yintercept=0, color="white", linetype="dashed") +
  theme_black() +
  xlim(c(0, 21.5)) + ylim(c(-1, 3.2)) +
  scale_colour_manual(values = rev(pal2), labels = c("Tropics n = 37", 
                                                     "Subtropics n = 64",
                                                     "Temperate n = 70", 
                                                     "Polar n = 2")) +
  scale_fill_manual(values = rev(pal2), labels = c("Tropics n = 37", 
                                                   "Subtropics n = 64",
                                                   "Temperate n = 70", 
                                                   "Polar n = 2")) 

#ggsave("Output/Figures/interaction_effects_trop_sub_temp", device = "png",
#    height = 9, width = 16, dpi = 400)




# Graph log interaction model ------
# make own df of lat values
mean_lat_chla_log <- meta_final %>%
  group_by(climate2) %>%
  summarise(mean = mean(abs_lat),
            min_chla_log = min(chla_log),
            max_chla_log = max(chla_log))

# what are the min and max chl a values for each zone
lat_values_log <- data.frame(chla_log = c(seq(from = -2.3280832, to = 0.4201407, length.out = 20),
                                          seq(from = -2.9682384, to = 2.3666108, length.out = 20),
                                          seq(from = -0.1570069, to = 3.0673489, length.out = 20),
                                          seq(from = -0.8074610, to = -0.8074610, length.out = 20)),
                             abs_lat = c(rep(17.45767, times = 20), 
                                         rep(27.79097, times = 20), 
                                         rep(47.36229, times = 20), 
                                         rep(74.21667, times = 20))) %>%
  mutate(`chla_log:abs_lat` = chla_log*abs_lat) %>%
  as.matrix()

# change the sequence for chl a to represent the chl a values we actually see at these latitudes (similar notation to how I did latitude with c() instead of rep())
# and change lat values for each bin to better reflect the avg lat in that bin

# use predict to tell it which values you want predictions for
lat_chla_predict_log <- as.data.frame(predict(object = model_multi_log,
                                              newmods = lat_values_log)) %>%
  cbind(lat_values_log) %>%
  mutate(climate = factor(c(rep("Tropics", times = 20), 
                            rep("Subtropics", times = 20), 
                            rep("Temperate", times = 20), 
                            rep("Polar", times = 20)),
                          levels = c("Tropics", "Subtropics", "Temperate", "Polar")))

# plot predicted values as line, jitter observed effect sizes
ggplot() +
  geom_line(data = lat_chla_predict_log,  aes(chla_log, pred, colour = climate)) +
  geom_point(data = meta_final, aes(chla_log, yi, colour = climate2)) +
  geom_ribbon(data = lat_chla_predict_log, aes(chla_log, pred, ymin = ci.lb, ymax = ci.ub, fill = climate), alpha = .1) +
  labs(x = "Log(chl a)", y = "Predicted Effect Size", colour = "Climate", fill = "Climate") +
  geom_hline(yintercept=0, color="white", linetype="dashed") +
  theme_black()


#ggsave("Output/Figures/interaction_effects_log", device = "png",
#     height = 9, width = 16, dpi = 400)

# Chl a x cat lat interaction model ------
model_multi2 <- rma.mv(yi = yi,
                       V = vi,
                       slab = author,
                       data = meta_final,
                       method = "REML",
                       mods = ~ chla_est * climate2,
                       random = ~ 1 | paper_num/indiv,
                       test = "t")

summary(model_multi2)

# Compare corals to endosymbionts? -----
corals <- meta_final %>% 
  filter(pp_taxa == "Co" |pp_taxa == "En")

model_coral <- rma.mv(yi = yi,
                      V = vi,
                      slab = author,
                      data = corals,
                      method = "REML",
                      mods = ~ pp_taxa,
                      random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                      test = "t")

summary(model_coral)
# Look at Test of Moderators to detemine whether there are significant differences between the subgroups
# p-value is 0.8314 so no, there's no diff between corals and endosymbionts
# The Model Results are printed within a meta-regression framework. This means that we cannot directly extract the estimates in order to obtain the pooled effect sizes within subgroups.

# How many papers in each bin?
coral_n <- corals %>%
  group_by(pp_taxa) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:2,
         pp_taxa = c("Coral", "Endosymbiont") )

# Pull coefficient estimates so I can graph them
intercept_coral <- coef(summary(model_coral))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes

# edit bin labels
coral_effect_sum <- coef(summary(model_coral)) %>%
  mutate(pp_taxa = factor(c("Coral", "Endosymbiont"),
                          levels = c("Coral", "Endosymbiont")),
         effect = ifelse(pp_taxa != "Coral", estimate + intercept_coral,
                         intercept_coral),
         ci.lb = ifelse(pp_taxa != "Coral", ci.lb + intercept_coral,
                        ci.lb),
         ci.ub = ifelse(pp_taxa != "Coral", ci.ub + intercept_coral,
                        ci.ub)) %>% # add intercept to other coefficients
  left_join(coral_n, by = "pp_taxa")

# GRAPH THIS?
coral_labs <- c("Coral", "Endosymbiont")

coral_n_labs <- c("n = 28", "n = 8")

ggplot(coral_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = pp_taxa)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed", alpha=.5) +
  theme_classic() +
  labs(x = "Effect size (ROM)", y = " ") +
  scale_y_continuous(breaks = 1:length(coral_labs),
                     labels = coral_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(coral_n_labs),
                                         labels = coral_n_labs)) +
  theme_black() +
  theme(legend.position = "none")
# Ok so yeah basically the coral and endosymbionts aren't showing different effect sizes, the coral is just a tighter estimate, probably bc it has a larger n

# Multi-level models with single mods or interactions -------

# Null model ------

model_null <- rma.mv(yi = yi,
                     V = vi,
                     slab = author,
                     data = meta_final,
                     method = "REML",
                     random = ~ 1 | paper_num/indiv, 
                     test = "t")

summary(model_null) # <.0001


# Latitude mulit level model ----

# We are going to use the absolute value of latitude, because I have no reason to expect North temperate to vary from South temperate
model_lat <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = meta_final,
                    method = "REML",
                    mods = ~ lat_stand,
                    random = ~ 1 | paper_num/indiv, 
                    test = "t")

summary(model_lat)
# meta_final --> p-val = 0.7132
# meta_umolL --> p-val = p-val = 0.1782

# Look at Test of Moderators to determine whether there are significant differences between the subgroups
# p-value is 0.6181 so NO, there is NO difference!
# The Model Results are printed within a meta-regression framework. This means that we cannot directly extract the estimates in order to obtain the pooled effect sizes within subgroups.

# Graph lat 

# Use predict() to predict effect size for each abs lat
# https://wviechtb.github.io/metafor/reference/predict.rma.html

# For random-effects models that were fitted with the rma.mv function, the model may actually include multiple 𝜏2 values (i.e., when the random argument includes an ‘~ inner | outer’ term and struct="HCS", struct="DIAG", struct="HAR", or struct="UN"). In that case, the function will provide prediction intervals for each level of the inner factor (since the prediction intervals differ depending on the 𝜏2 value)

# don't have to worry about tau2.levels and gamma2.levels arguments because I'm not using newmods, I'm just predicting the effects at the lat values from my input data

lat_predict <- as.data.frame(predict(model_lat,
                                     addx=TRUE)) %>%
  tibble::rownames_to_column("author") %>%
  rename(abs_lat = X.abs_lat) 

# pred = predicted value(s).
# se = corresponding standard error(s).
#ci.lb = lower bound of the confidence interval(s).
#ci.ub = upper bound of the confidence interval(s).
#pi.lb = lower bound of the prediction interval(s) (only for random/mixed-effects models).
#pi.ub = upper bound of the prediction interval(s) (only for random/mixed-effects models).

# Use ggeffects to graph these predicted values
# https://strengejacke.github.io/ggeffects/index.html
# geom line for predicted and geom ribbon for CI intervals

# Create a palette
# how many unique papers
n <- length(unique(meta_final$paper_num))

#combine the qualitative colour brewer palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# set seed because otherwise I'll get different colours each time
set.seed(1235)

# plot predicted values as line, jitter observed effect sizes
ggplot() +
  geom_line(data = lat_predict, aes(abs_lat, pred), colour = "white", size = 1.5) +
  geom_jitter(data = meta_final, aes(abs_lat, yi, colour = author)) +
  geom_ribbon(data = lat_predict, aes(abs_lat, pred, ymin = ci.lb, ymax = ci.ub), alpha = .1, fill = "white") +
  theme_classic() +
  labs(x = "Absolute latitude", y = "Effect Size") +
  scale_color_manual(values = sample(col_vector, n)) +
  geom_hline(yintercept=0, color="white", linetype="dashed") +
  theme_black() +
  theme(legend.position="none")

ggsave("Output/Figures/abs_lat_effect", device = "png",
       height = 9, width = 16, dpi = 400)


# Background nutrients multilevel model-----
model_back_nutrients <- rma.mv(yi = yi,
                               V = vi,
                               slab = author,
                               data = meta_final,
                               method = "REML",
                               mods = ~ background_nut,
                               random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                               test = "t")

summary(model_back_nutrients)
# meta_final --> p-value = 0.2242
# meta_umolL --> p-value = 0.0162

# Graph background nutrients 
nut_predict <- as.data.frame(predict(model_back_nutrients,
                                     addx=TRUE)) %>%
  tibble::rownames_to_column("author") %>%
  rename(background_nut = X.background_nut) 

# set seed because otherwise I'll get different colours each time
set.seed(1235)

# plot predicted values as line, jitter observed effect sizes
ggplot() +
  geom_line(data = nut_predict, aes(background_nut, pred), colour = "white", size = 1.5) +
  geom_jitter(data = meta_final, aes(background_nut, yi, colour = climate2)) +
  geom_ribbon(data = nut_predict, aes(background_nut, pred, ymin = ci.lb, ymax = ci.ub), fill = "white", alpha = .1) +
  theme_classic() +
  labs(x = "Background nutrients (umol/L)", y = "Effect Size", colour = "Climate Zone") +
  geom_hline(yintercept=0, color="white", linetype="dashed") +
  theme_black() +
  scale_colour_manual(values = rev(pal2))

ggsave("Output/Figures/background_nut_effect", device = "png",
       height = 9, width = 16, dpi = 400)

# Added nitrogen multilevel model ------
# Only looking at added nitrogen where they measured it in umol/L
meta_umolL <- meta_final %>%
  filter(prov_units == "umol/L") %>%
  filter(added_N < 10) # Filter out the super high ammonium estimates that are dragging the line down?

model_added_N <- rma.mv(yi = yi,
                        V = vi,
                        slab = author,
                        data = meta_umolL,
                        method = "REML",
                        mods = ~ added_N,
                        random = ~ 1 | paper_num/indiv, 
                        test = "t")

summary(model_added_N)
# p-val = 0.1853

model_added_nh4 <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_umolL,
                          method = "REML",
                          mods = ~ added_nh4,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_added_nh4)
# p-val = 0.5276

# Graph added nitrogen 
N_predict <- as.data.frame(predict(model_added_N,
                                   addx=TRUE)) %>%
  tibble::rownames_to_column("author") %>%
  rename(added_N = X.added_N) 

# set seed because otherwise I'll get different colours each time
set.seed(1235)

# plot predicted values as line, jitter observed effect sizes
ggplot() +
  geom_line(data = N_predict, aes(added_N, pred), colour = "white", size = 1.5) +
  geom_jitter(data = meta_umolL, aes(added_N, yi, colour = climate2)) +
  geom_ribbon(data = N_predict, aes(added_N, pred, ymin = ci.lb, ymax = ci.ub), fill = "white", alpha = .1) +
  theme_classic() +
  labs(x = "Added nitrogen (umol/L)", y = "Effect Size", colour = "Climate Zone") +
  geom_hline(yintercept=0, color="white", linetype="dashed") +
  theme_black() +
  scale_colour_manual(values = rev(pal2))

#ggsave("Output/Figures/added_N_effect", device = "png",
# height = 9, width = 16, dpi = 400)

# This shows that for the most part, added N was around 1-10
# Should I just cut the outliers that were above 10..? It's just that one mussel vs oyster vs salmon vs bass paper, the three high points are three of the four experiments
meta_big_pee <- meta_umolL %>%
  filter(added_N > 10)

# Temperature multilevel model-----
# Temperature and abs_lat are correlated 5.0679
model_temp <- rma.mv(yi = yi,
                     V = vi,
                     slab = author,
                     data = meta_final,
                     method = "REML",
                     mods = ~ temp,
                     random = ~ 1 | paper_num/indiv,
                     test = "t")

summary(model_temp)
# meta_final --> p-val = 0.5513
# meta_umolL --> p-val = 0.0995

# Graph temperature
temp_predict <- as.data.frame(predict(model_temp,
                                      addx=TRUE)) %>%
  tibble::rownames_to_column("author") %>%
  rename(temp = X.temp) 

# set seed because otherwise I'll get different colours each time
set.seed(1235)

# plot predicted values as line, jitter observed effect sizes
ggplot() +
  geom_line(data = temp_predict, aes(temp, pred), colour = "white", size = 1.5) +
  geom_jitter(data = meta_final, aes(temp, yi, colour = climate2)) +
  geom_ribbon(data = temp_predict, aes(temp, pred, ymin = ci.lb, ymax = ci.ub), fill = "white", alpha = .1) +
  theme_classic() +
  labs(x = "Temperature (C)", y = "Effect Size", colour = "Climate Zone") +
  geom_hline(yintercept=0, color="white", linetype="dashed") +
  theme_black() +
  scale_colour_manual(values = rev(pal2))

#ggsave("Output/Figures/temperature_effect", device = "png",
#     height = 9, width = 16, dpi = 400)


# Experimental design multilevel model-----
model_exp_design <- rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = meta_final,
                           method = "REML",
                           mods = ~ natural,
                           random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                           test = "t")

summary(model_exp_design)
# meta_final --> p-val = 0.1796
# meta_umolL --> p-val = 0.2430

# How many papers in each bin?
nat_n <- meta_final %>%
  group_by(natural) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:5,
         nat_bin = factor(c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm"),
                          levels = c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")))

# Pull coefficient estimates so I can graph them
intercept_nat <- coef(summary(model_exp_design))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes

# When you add intercept to the other estimates, the p-value asks whether the CI’s intercept the intercept, not whether they intercept 0

# Can’t just add intercept to CI, works for the estimates but not the CI
# Try using emmeans for this

# edit bin labels
nat_effect_sum <- coef(summary(model_exp_design)) %>%
  mutate(nat_bin = factor(c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm"),
                          levels = c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")),
         effect = ifelse(nat_bin != "Field: unstocked", estimate + intercept_nat,
                         intercept_nat),
         ci.lb = ifelse(nat_bin != "Field: unstocked", ci.lb + intercept_nat,
                        ci.lb),
         ci.ub = ifelse(nat_bin != "Field: unstocked", ci.ub + intercept_nat,
                        ci.ub)) %>% # add intercept to other coefficients
  left_join(nat_n, by = "nat_bin")

# GRAPH THIS?
nat_labs <- c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")

nat_n_labs <- c("n = 74", "n = 14", "n = 21", "n = 9", "n = 55")

ggplot(nat_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = nat_bin)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size", y = " ") +
  scale_y_continuous(breaks = 1:length(nat_labs),
                     labels = nat_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(nat_n_labs),
                                         labels = nat_n_labs)) +
  scale_colour_manual(values = rev(pal_star)) +
  theme_black()+ 
  theme(legend.position = "none")

#ggsave("Output/Figures/effects_nat", device = "png",
#     height = 9, width = 8, dpi = 400)

# Consumer activity level multilevel model -----
model_act <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = meta_final,
                    method = "REML",
                    mods = ~ c_activity,
                    random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                    test = "t")

summary(model_act)
# meta_final --> p-value = 0.2224 
# meta_umolL --> p-value = 0.0982 

# Graph consumer level 
# How many papers in each bin?
act_n <- meta_final %>%
  group_by(c_activity) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:3)

# Pull coefficient estimates so I can graph them
intercept_act <- coef(summary(model_act))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes

act_effect_sum <- coef(summary(model_act)) %>%
  mutate(c_activity = factor(c("Vertebrate", "Mobile Invert", "Sessile Invert"),
                             levels = c("Vertebrate", "Mobile Invert", "Sessile Invert")),
         effect = ifelse(c_activity != "Vertebrate", estimate + intercept_act,
                         intercept_act),
         ci.lb = ifelse(c_activity != "Vertebrate", ci.lb + intercept_act,
                        ci.lb),
         ci.ub = ifelse(c_activity != "Vertebrate", ci.ub + intercept_act,
                        ci.ub)) %>% # add intercept to other coefficients
  left_join(act_n, by = "c_activity")

# GRAPH THIS?
act_labs <- c("Vertebrate", "Mobile Invert", "Sessile Invert")

act_n_labs <- c("n = 55", "n = 20", "n = 98")

ggplot(act_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = c_activity)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size", y = " ") +
  scale_y_continuous(breaks = 1:length(act_labs),
                     labels = act_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(act_n_labs),
                                         labels = act_n_labs)) +
  scale_colour_manual(values = rev(pal_star)) +
  theme_black()+ 
  theme(legend.position = "none")

#ggsave("Output/Figures/effects_act", device = "png",
#       height = 9, width = 8, dpi = 400)


# I feel like at some point I will have to address the differences between the tropical taxons and the temperate taxons

# Does activity level predict provisioned nutrients?

meta_estimated <- meta_final %>%
  filter(nutrient_estimate == "measured")

model_activity_lmer <- lmer(added_nh4 ~ c_activity + (1|paper_num), data = meta_estimated)
summary(model_activity_lmer)
# Not really

# Primary producer multilevel model -----
# Groups are Algae (FA + CA), Coral (En + Co + anemone), Phytoplankton (P), and Plants (S + mangrove)
model_pp <- rma.mv(yi = yi,
                   V = vi,
                   slab = author,
                   data = meta_final,
                   method = "REML",
                   mods = ~ pp_group ,
                   random = ~ 1 | paper_num/indiv,
                   test = "t")

summary(model_pp) 
# meta_final --> p-val = 0.0278
# meta_umolL --> p-val = 0.0173

# OK so producer ID matters, algae show larger effects vs the rest
# Ok phytoplankton are the best distributed let's just look at them for a sec
algae <- meta_final %>%
  filter(pp_group == "Algae")
# I just looked at the phytoplankton and the algae
# Both had the same trend

# OKKKK does algae effect vary by lat
p_lat <- rma.mv(yi = yi,
                V = vi,
                slab = author,
                data = algae,
                method = "REML",
                mods = ~ abs_lat,
                random = ~ 1 | paper_num/indiv,
                test = "t")

summary(p_lat) # NO, p-val = 0.2787

# OKKKK does algae effect vary by back nut
p_back <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = algae,
                 method = "REML",
                 mods = ~ background_nut,
                 random = ~ 1 | paper_num/indiv,
                 test = "t")

summary(p_back) # YES, p-val = 0.0738
# Larger effect of pee on algae with higher levels of background nutrients

# OKKKK does algae effect vary by added nut or added nh4
# Use meta_umolL because that only has the measured nitrogen added in umol/L

algae_umol <- meta_umolL %>%
  filter(pp_group == "Algae")
#sample size is only 10 we can't do this

p_added_N <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = algae_umol,
                    method = "REML",
                    mods = ~ added_N,
                    random = ~ 1 | paper_num/indiv,
                    test = "t")

summary(p_added_N) # NO, p-val = 0.7897

# OKKKK does algae effect vary by temperature
p_temp <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = algae,
                 method = "REML",
                 mods = ~ temp,
                 random = ~ 1 | paper_num/indiv,
                 test = "t")

summary(p_temp) # NO, p-val = 0.1676

# OK SO 

# Graph primary producer
# How many papers in each bin?
pp_n <- meta_final %>%
  group_by(pp_group) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:4)

# Pull coefficient estimates so I can graph them
intercept_pp <- coef(summary(model_pp))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes

pp_effect_sum <- coef(summary(model_pp)) %>%
  mutate(pp_group = factor(c("Algae", "Coral", "Phytoplankton", "Plants")), 
         effect = ifelse(pp_group != "Algae", estimate + intercept_pp,
                         intercept_pp),
         ci.lb = ifelse(pp_group != "Algae", ci.lb + intercept_pp,
                        ci.lb),
         ci.ub = ifelse(pp_group != "Algae", ci.ub + intercept_pp,
                        ci.ub)) %>% # add intercept to other coefficients
  left_join(pp_n, by = "pp_group")

# GRAPH THIS?
pp_labs <- c("Algae", "Coral", "Phytoplankton", "Plants")

pp_n_labs <- c("n = 26", "n = 47", "n = 53", "n = 47")

ggplot(pp_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = pp_group)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept = intercept_pp, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size", y = " ") +
  scale_y_continuous(breaks = 1:length(pp_labs),
                     labels = pp_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(pp_n_labs),
                                         labels = pp_n_labs)) +
  scale_colour_manual(values = rev(pal_star)) +
  theme_black()+ 
  theme(legend.position = "none")

#ggsave("Output/Figures/effects_pp", device = "png",
#       height = 9, width = 16, dpi = 400)


# Interaction latitude x background nutrients -----
model_lat_back <- rma.mv(yi = yi,
                         V = vi,
                         slab = author,
                         data = meta_final,
                         method = "REML",
                         mods = ~ background_nut * abs_lat,
                         random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                         test = "t")

summary(model_lat_back)
vif(model_lat_back)
# meta_final --> p-val = 0.0480
# meta_umolL --> p-val = 0.1055

# So overall the test of mods is almost signif (p-val = 0.1079)
# but each mod or interaction isn't signif
# when you run the model without any interactions background nut is positive and signif and lat is negative and signif
# so as you increase background nutrients the effect increases
# so as you increase lat the effect decreases

# check for colinearity
vif(model_lat_back) # run this on the model without an interaction
# value of 1 means they're not correlated with any other variables
# higher value = greater correlation
# values of more than 4 or 5 are regarded as moderate - high
# values of 10+ are very high

# OK so there's low colinearity (2.2405)

# What is the relationship between the two
lm_nut_lat <- lmer(background_nut ~ abs_lat + (1|paper_num), data = meta_final)
summary(lm_nut_lat)
visreg(lm_nut_lat)
# there's a link but the R2 is quite low (0.11)

# Graphing interaction between nut and lat 
# make own df of lat values
mean_lat_nut <- meta_final %>%
  group_by(climate2) %>%
  summarise(mean = mean(abs_lat),
            min_nut = min(background_nut, na.rm = TRUE),
            max_nut = max(background_nut, na.rm = TRUE))

# what are the min and max chl a values for each zome


lat_back_values <- data.frame(background_nut = c(seq(from = 0.00000117, to = 1.58000, length.out = 20),
                                                 seq(from = 0.06521739, to = 9.96201, length.out = 20),
                                                 seq(from = 1, to = 23, length.out = 20),
                                                 seq(from = 10.61666667, to = 10.61666667, length.out = 20)),
                              abs_lat = c(rep(17.42045, times = 20), 
                                          rep(26.54635, times = 20), 
                                          rep(48.00235, times = 20), 
                                          rep(74.21667, times = 20))) %>%
  mutate(`background_nut:abs_lat` = background_nut*abs_lat) %>%
  as.matrix()

# change the sequence for chl a to represent the chl a values we actually see at these latitudes (similar notation to how I did latitude with c() instead of rep())
# and change lat values for each bin to better reflect the avg lat in that bin

# use predict to tell it which values you want predictions for
lat_back_predict <- as.data.frame(predict(object = model_lat_back,
                                          newmods = lat_back_values)) %>%
  cbind(lat_back_values) %>%
  mutate(climate = factor(c(rep("Tropics", times = 20), 
                            rep("Subtropics", times = 20), 
                            rep("Temperate", times = 20), 
                            rep("Polar", times = 20)),
                          levels = c("Tropics", "Subtropics", "Temperate", "Polar")))

# determine bin sizes
bin_size <- meta_final %>%
  drop_na(background_nut) %>%
  filter(climate2 == "Temperate")

# plot predicted values as line, jitter observed effect sizes
ggplot() +
  geom_line(data = lat_back_predict,  aes(background_nut, pred, colour = climate), size = 1.7) +
  geom_point(data = meta_final, aes(background_nut, yi, colour = climate2), size = 2.5) +
  geom_ribbon(data = lat_back_predict, aes(background_nut, pred, ymin = ci.lb, ymax = ci.ub, fill = climate), alpha = .2) +
  labs(x = "Background nutrients (umol/L)", y = "Effect Size", colour = "Climate", fill = "Climate") +
  geom_hline(yintercept=0, color="white", linetype="dashed") +
  theme_black() +
  xlim(c(0, 21.5)) + ylim(c(-1, 3.2)) +
  scale_colour_manual(values = rev(pal2), labels = c("Tropics n = 37", 
                                                     "Subtropics n = 64",
                                                     "Temperate n = 69", 
                                                     "Polar n = 2")) +
  scale_fill_manual(values = rev(pal2), labels = c("Tropics n = 37", 
                                                   "Subtropics n = 64",
                                                   "Temperate n = 69", 
                                                   "Polar n = 2")) 

#ggsave("Output/Figures/interaction_nut_lat_effects", device = "png",
#     height = 9, width = 16, dpi = 400)

# Interaction added nutrients x background nutrients ----
# Use meta_umolL because that only has the measured nitrogen added in umol/L
model_added_background <- rma.mv(yi = yi,
                                 V = vi,
                                 slab = author,
                                 data = meta_umolL,
                                 method = "REML",
                                 mods = ~ added_N * background_nut,
                                 random = ~ 1 | paper_num/indiv, 
                                 test = "t")

summary(model_added_background)
# p value = 0.1075
vif(model_added_background) # 1.1285 

# Interaction temperature x added nutrients ------
# Use meta_umolL because that only has the measured nitrogen added in umol/L
model_temp_added <- rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = meta_umolL,
                           method = "REML",
                           mods = ~ temp * added_N,
                           random = ~ 1 | paper_num/indiv, 
                           test = "t")

summary(model_temp_added)
# p-val = 0.3175

vif(model_temp_added) #  1.0782 

# Interaction temperature x background nutrients ------
model_temp_back <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_final,
                          method = "REML",
                          mods = ~ temp * background_nut,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_temp_back)
# meta_final --> p-val = 0.4587
# meta_umolL --> p-val = 0.1126

vif(model_temp_back) # 1.9817


# AIC model comparisons -------
# Careful here bc each model has different data
# OK I went through and ran each model on the same data (meta_umolL)
# And recorded the differences in p-values
# And the model with added_N and background_nutrients is the best


AIC_values <- AIC(model_null, model_lat, model_back_nutrients, model_added_N, model_temp, model_exp_design, model_act, model_pp, model_lat_back, model_added_background, model_temp_added, model_temp_back)
# Model with added and background is the best?????????

BIC_values <- BIC(model_null, model_lat, model_back_nutrients, model_added_N, model_temp, model_exp_design, model_act, model_pp, model_lat_back, model_added_background, model_temp_added, model_temp_back)
# Model with added N and model with added N and background is best here

# OK SHIT SHIT SHIT now when I use the meta_umolL dataset for all the models the AIC is all jumbled
# And lots of things are significant that were not before
# OH GOD



# Influential studies with single mod models --------
meta_drop <- meta_final %>%
  filter(paper_num != 17) 

meta_umolL_drop <- meta_umolL %>%
  filter(paper_num != 17) %>%
  filter(paper_num != 23)

# Summary of drop one at a time changes:
# consumer activity = 23 makes it signif (0.0183 vs 0.1677)
# lat x background = 17 makes it signif (0.0480 vs 0.1079) 
# added_N x background = 17 makes it almost signif  (0.1075 vs 0.3726) 
# PP = 17 makes it slightly more signif (0.0278 vs 0.0601)
# Experimental design = 29 makes it not signif but closer..? (0.0982 vs 0.1998)

# When you drop both 17 and 23 these mods become signif:
# background_nutrients (0.0462 versus 0.3141)
# consumer activity level (0.0268 versus 0.1677)
# added_n x background (0.1075 versus 0.3726)

# OK HERE'S THE DECISION. I DON'T THINK CONSUMER ACTIVITY LEVEL ACTUALLY MATTERS. SO I'M NOT DROPPING 23 JUST BECAUSE IT AFFECTS CONSUMER ACTIVITY LEVEL. SO JUST GONNA DROP 17 (CAVES WITH SPONGES)

# Rerun lat ----
re_model_lat <- rma.mv(yi = yi,
                       V = vi,
                       slab = author,
                       data = meta_drop,
                       method = "REML",
                       mods = ~ abs_lat,
                       random = ~ 1 | paper_num/indiv, 
                       test = "t")

summary(re_model_lat)

# Removed # 17, p-value is 0.7132 versus 0.6181
# Removed # 23, p-value is 0.6927 versus 0.6181
# Removed # 29, p-value is 0.3172 versus 0.6181
# Removed # 17 and 23, p-value is 0.5903 versus 0.6181

# Rerun background nutrients-----
re_model_back_nutrients <- rma.mv(yi = yi,
                                  V = vi,
                                  slab = author,
                                  data = meta_drop,
                                  method = "REML",
                                  mods = ~ background_nut,
                                  random = ~ 1 | paper_num/indiv,
                                  test = "t")

summary(re_model_back_nutrients)

# Remove # 17 p-value is 0.2242 versus p-val =  0.3141
# Remove # 23 p-value is 0.1181 versus p-val =  0.3141
# Remove # 29 p-value is 0.7017 versus p-val =  0.3141
# Removed # 17 and 23, p-value is 0.0462 versus 0.3141

# Rerun added nitrogen ------
# Only looking at added nitrogen where they measured it in umol/L
re_model_added_N <- rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = meta_umolL_drop,
                           method = "REML",
                           mods = ~ added_N,
                           random = ~ 1 | paper_num/indiv, 
                           test = "t")

summary(re_model_added_N)

# Remove # 17 p-value is 0.1853 versus p-val =  0.1976
# Remove # 23 p-value is 0.1976 versus p-val =  0.1976
# Remove # 29 p-value is 0.2290 versus p-val =  0.1976
# Removed # 17 and 23, p-value is 0.1853 versus 0.1976

# Rerun temperature multilevel model-----
re_model_temp <- rma.mv(yi = yi,
                        V = vi,
                        slab = author,
                        data = meta_drop,
                        method = "REML",
                        mods = ~ temp,
                        random = ~ 1 | paper_num/indiv,
                        test = "t")

summary(re_model_temp)

# Remove # 17 p-value is 0.5513 versus p-val =  0.6678
# Remove # 23 p-value is 0.2581 versus p-val =  0.6678
# Remove # 29 p-value is 0.9791 versus p-val =  0.6678
# Removed # 17 and 23, p-value is 0.1525 versus 0.6678

# Rerun experimental design -----
re_model_exp_design <- rma.mv(yi = yi,
                              V = vi,
                              slab = author,
                              data = meta_drop,
                              method = "REML",
                              mods = ~ natural,
                              random = ~ 1 | paper_num/indiv,
                              test = "t")

summary(re_model_exp_design)

# Remove # 17 p value is 0.1796 vs 0.1998 
# Remove # 23 p value is 0.2493 vs 0.1998 
# Remove # 29 p value is 0.0982 vs 0.1998 
# Removed # 17 and 23, p-value is 0.3051 versus 0.1998

# Rerun consumer activity level -----
re_model_act <- rma.mv(yi = yi,
                       V = vi,
                       slab = author,
                       data = meta_drop,
                       method = "REML",
                       mods = ~ c_activity,
                       random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                       test = "t")

summary(re_model_act)

# Removed # 17 p-value is 0.2224 vs 0.1677 
# Removed # 23 p-value is 0.0183 vs 0.1677 
# Removed # 29 p-value is 0.2244 vs 0.1677 
# Removed # 17 and 23, p-value is 0.0268 versus 0.1677

# Rerun primary producer -----
re_model_pp <- rma.mv(yi = yi,
                      V = vi,
                      slab = author,
                      data = meta_drop,
                      method = "REML",
                      mods = ~ pp_group ,
                      random = ~ 1 | paper_num/indiv,
                      test = "t")

summary(re_model_pp) 

# Removed # 17 p-value is 0.0278 vs 0.0601 
# Removed # 23 p-value is 0.0627 vs 0.0601 
# Removed # 29 p-value is 0.0093 vs 0.0601 
# Removed # 17 and 23, p-value is 0.0244 versus 0.0601


# Rerun latitude x background nutrients -----
re_model_lat_back <- rma.mv(yi = yi,
                            V = vi,
                            slab = author,
                            data = meta_drop,
                            method = "REML",
                            mods = ~ background_nut * abs_lat,
                            random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                            test = "t")

summary(re_model_lat_back)

# Removed # 17 p-value is 0.0480 vs 0.1079 
# Removed # 23 p-value is 0.2204 vs 0.1079 
# Removed # 29 p-value is 0.2186 vs 0.1079 
# Removed # 17 and 23, p-value is 0.0892 versus 0.1079

# Rerun added nutrients x background nutrients ----
# Use meta_umolL because that only has the measured nitrogen added in umol/L
re_model_added_background <- rma.mv(yi = yi,
                                    V = vi,
                                    slab = author,
                                    data = meta_umolL_drop,
                                    method = "REML",
                                    mods = ~ added_N*background_nut,
                                    random = ~ 1 | paper_num/indiv, 
                                    test = "t")

summary(re_model_added_background)

# Removed # 17 p-value is 0.1075 vs 0.3726 
# Removed # 23 p-value is 0.3726 vs 0.3726 
# Removed # 29 p-value is 0.2793 vs 0.3726 
# Removed # 17 and 23, p-value is 0.1075 versus 0.3726


# Rerun temperature x added nutrients ------
# Use meta_umolL because that only has the measured nitrogen added in umol/L
re_model_temp_added <- rma.mv(yi = yi,
                              V = vi,
                              slab = author,
                              data = meta_umolL_drop,
                              method = "REML",
                              mods = ~ temp * added_N,
                              random = ~ 1 | paper_num/indiv, 
                              test = "t")

summary(re_model_temp_added)

# Removed # 17 p-value is 0.3175 vs 0.5870 
# Removed # 23 p-value is 0.5870 vs 0.5870 
# Removed # 29 p-value is 0.6790 vs 0.5870 
# Removed # 17 and 23, p-value is 0.3175 versus 0.5870


# Rerun temperature x background nutrients ------
re_model_temp_back <- rma.mv(yi = yi,
                             V = vi,
                             slab = author,
                             data = meta_drop,
                             method = "REML",
                             mods = ~ temp * background_nut,
                             random = ~ 1 | paper_num/indiv, 
                             test = "t")

summary(re_model_temp_back)

# Removed # 17 p-value is 0.4587 vs 0.4673 
# Removed # 23 p-value is 0.2985 vs 0.4673 
# Removed # 29 p-value is 0.8587 vs 0.4673 
# Removed # 17 and 23, p-value is 0.2063 versus 0.4673

AIC(re_model_act, re_model_added_background, re_model_added_N, re_model_back_nutrients, re_model_exp_design, re_model_lat, re_model_lat_back, re_model_pp, re_model_temp, re_model_temp_added, re_model_temp_back)
# Added x background still the best model

temp <- meta_final %>%
  filter(abs_lat > 30) %>%
  filter(obs_num == 1)




# old stats -----
# Full versus reduced model comparisons in terms of fit statistics and likelihoods can be obtained with anova.rma
# Permutation tests for the model coefficient(s) can be obtained with permutest.rma.uni. Tests and confidence intervals based on (cluster) robust methods can be obtained with robust.rma.uni.
# Tests for funnel plot asymmetry (which may be indicative of publication bias) can be obtained with ranktest.rma and regtest.rma. For models without moderators, the trimfill.rma.uni method can be used to carry out a trim and fill analysis and hc.rma.uni provides a random-effects model analysis that is more robust to publication bias (based on the method by Henmi & Copas, 2010).

# Try using "Doing Meta-Analysis in R"
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/metareg.html

# If you don't want to compare the models directly using the anova function, use the "REML" (restricted maximum likelihood) 
model_control <- rma(yi = yi,
                     vi = vi,
                     data = meta_final,
                     method = "ML",
                     test = "knha")

model_lat <- rma(yi = yi,
                 vi = vi,
                 data = meta_final,
                 method = "ML",
                 mods = ~ lat_dec,
                 test = "knha")

# is the model with lat better than the null model? NO
anova(model_control, model_lat)

model_chla <- rma(yi = yi,
                  vi = vi,
                  data = meta_final,
                  method = "ML",
                  mods = ~ chla_bin,
                  test = "knha")

# how about the model with chl a? is it better than the null model? YES???
anova(model_control, model_chla)

# Look at effect sizes????
ggplot(meta, aes(x = yi, y = lat_dec,  colour = climate_zone)) +
  geom_point(size = 3) +
  scale_colour_manual(values = rev(pal), labels = climate_labs, guide = guide_legend(reverse = TRUE)) +
  labs( y = "Lattitude", x = "Effect size", colour = "Climate Zone") + 
  geom_vline(xintercept = 0, linetype="dotted", colour = "white", size = 1.6) +
  theme_black() +
  ylim(c(-75, 75))

#ggsave("Output/Figures/response_lat", device = "png",
#       height = 9, width = 16, dpi = 400)



ggplot(meta, aes(x = yi, y = lat_dec,  colour = climate_zone)) +
  geom_point(size = 3) +
  scale_colour_manual(values = rev(pal), guide = "none") +
  labs( y = "Lattitude", x = "Effect size", colour = "Climate Zone") + 
  geom_vline(xintercept = 0, linetype="dotted", colour = "white", size = 1.6) +
  theme_black()

#ggsave("Output/Figures/response_lat_noleg", device = "png",
#     height = 9, width = 16, dpi = 400)

# Plotting binned chla ---------
# How many papers in each bin?
chla_n <- meta_final %>%
  group_by(chla_bin) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:5,
         chla_bin = factor(c("Very low", "Low", "Medium", "High", "Very high")))

# Pull coefficient estimates so I can graph them
intercept_chla <- coef(summary(model_chla_paper_bin))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes
chla_effect_sum <- coef(summary(model_chla_paper_bin)) %>%
  mutate(chla_bin = factor(c("Very low", "Low", "Medium", "High", "Very high"),
                           levels = c("Very low", "Low", "Medium", "High", "Very high")),
         effect = ifelse(chla_bin != "Very low", estimate + intercept_chla,
                         intercept_chla),
         ci.lb = ifelse(chla_bin != "Very low", ci.lb + intercept_chla,
                        ci.lb),
         ci.ub = ifelse(chla_bin != "Very low", ci.ub + intercept_chla,
                        ci.ub)) %>%# add intercept to other coefficients
  left_join(chla_n, by = "chla_bin")

# GRAPH THIS?
chla_bin_labs <- c("Very low", "Low", "Medium", "High", "Very high")
n_labs_chla <- c("n = 109", "n = 29", "n = 7", "n = 6", "n = 18")

ggplot(chla_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = chla_bin)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept = 0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size (ROM)", y = "Chlorophyll a") +
  scale_y_continuous(breaks = 1:length(chla_bin_labs),
                     labels = chla_bin_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(n_labs_chla),
                                         labels = n_labs_chla)) +
  scale_color_brewer(palette = 'Greens') +
  theme_black()+ 
  theme(legend.position = "none")

#ggsave("Output/Figures/effects_chla", device = "png",
#      height = 9, width = 16, dpi = 400)


# Old graphing climate zones -------

# # Pull coefficient estimates so I can graph them
intercept_lat <- coef(summary(model_lat))[1,1] # extract the intercept

# How many papers in each bin?
lat_n <- meta_final %>%
  group_by(climate_zone) %>%
  summarise(n = n()) %>%
  rename(climate = climate_zone) %>%
  mutate(dummy = 1:6)

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes
lat_effect_sum <- coef(summary(model_lat)) %>%
  mutate(climate = factor(c("ST", "SST", "T", "NST", "NT", "NP"),
                          levels = c("ST", "SST", "T", "NST", "NT", "NP")),
         effect = ifelse(climate != "ST", estimate + intercept_lat,
                         intercept_lat),
         ci.lb = ifelse(climate != "ST", ci.lb + intercept_lat,
                        ci.lb),
         ci.ub = ifelse(climate != "ST", ci.ub + intercept_lat,
                        ci.ub)) %>%# add intercept to other coefficients
  left_join(lat_n, by = "climate")

# GRAPH THIS
n_labs_lat <- c("n = 5", "n = 10", "n = 28", "n = 50", "n = 79", "n = 2")

ggplot(lat_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = climate)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size (ROM)", y = "Climate zone") +
  scale_y_continuous(breaks = 1:length(climate_labs),
                     labels = climate_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(n_labs_lat),
                                         labels = n_labs_lat)) +
  scale_colour_manual(values = rev(pal))  +
  theme_black()+ 
  theme(legend.position = "none")

#ggsave("Output/Figures/effects_lat", device = "png",
#      height = 9, width = 16, dpi = 400)

# # Try using the 4 climate zones and graph ----
model_climate_paper <- rma.mv(yi = yi,
                              V = vi,
                              slab = author,
                              data = meta_final,
                              method = "REML",
                              mods = ~ climate2,
                              random = ~ 1 | paper_num/indiv, 
                              test = "t")

summary(model_climate_paper)
# Look at Test of Moderators to detemine whether there are significant differences between the subgroups
# p-value is 0.6538 so NO, there is NO difference!
# The Model Results are printed within a meta-regression framework. This means that we cannot directly extract the estimates in order to obtain the pooled effect sizes within subgroups.

# # Pull coefficient estimates so I can graph them
intercept_cli <- coef(summary(model_climate_paper))[1,1] # extract the intercept

# How many papers in each bin?
cli_n <- meta_final %>%
  group_by(climate2) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:4)

# Now pull the rest of the coefficients, add a column for climate, and adjust the estimates to reflect the estimated effect sizes
cli_effect_sum <- coef(summary(model_climate_paper)) %>%
  mutate(climate2 = factor(c("Tropics", "Subtropics", "Temperate", "Polar"),
                           levels = c("Tropics", "Subtropics", "Temperate", "Polar")),
         effect = ifelse(climate2 != "Tropics", estimate + intercept_cli,
                         intercept_cli),
         ci.lb = ifelse(climate2 != "Tropics", ci.lb + intercept_cli,
                        ci.lb),
         ci.ub = ifelse(climate2 != "Tropics", ci.ub + intercept_cli,
                        ci.ub)) %>% # add intercept to other coefficients
  left_join(cli_n, by = "climate2")


# GRAPH THIS
n_labs_cli <- c("n = 28", "n = 68", "n = 36", "n = 2")
cli_labs <- c("Tropics", "Subtropics", "Temperate", "Polar")

ggplot(cli_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = climate2)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size (ROM)", y = "Climate zone") +
  scale_y_continuous(breaks = 1:length(cli_labs),
                     labels = cli_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(n_labs_cli),
                                         labels = n_labs_cli)) +
  scale_colour_manual(values = rev(pal))  +
  theme_black()+ 
  theme(legend.position = "none")

#ggsave("Output/Figures/effects_climate", device = "png",
#      height = 9, width = 16, dpi = 400)


# Old way to get shapefiles ------
# Load chl a
#setwd("~/Downloads")
raster_summer2004 <- raster("Data/Rasters/A20041732004264.L3m_SNSU_CHL_chlor_a_4km.nc")
#setwd("~/Documents/MSc/Ch2-Meta-analysis")

# make a df for this one season / year
meta_summer2004 <- meta %>%
  filter(season == "summer" & year_of_study == "2004") %>%
  dplyr::select(indiv, lat_dec, long_dec)

write_csv(meta_summer2004, path =  "Output/Data/season_year/meta_summer2004.csv")

# Exported the meta_summer2004 df to csv, loaded it into QGIS, added a layer of points ('create points layer from table' tool), then exported it as a shapefile

# Open QGIS
# Layer > Add Layer > Add Deliminated Text Layer > select file
# File format = CSV
# Geometry Definition: Select Point coordinates, set x = long and y = lat
# Layer Settings: select Use spatial index
# Click Add
# Right click layer > Export > Save feature as > ESRI Shapefile
# Click dots next to file name to set path (put in shapefile folder)
# Click ok

sf_summer2004 <- readOGR(dsn = path.expand("~/Documents/MSc/Ch2-Meta-analysis/Shapefiles"), layer = "summer2004")

chla_summer2004 <- raster::extract(raster_summer2004, sf_summer2004[,c("lat_dec", "long_dec")], df = TRUE)

summer2004 <- cbind(meta_summer2004, chla_summer2004) %>%
  select(-ID)


# Try using tutorial to do meta stats? -----
#https://www.selfmindsociety.com/post/a-forest-plot-in-ggplot2
meta_model <- rma(yi, vi, data = meta )

meta_model2 <- rma(yi, vi, mods = ~ climate_zone, data = meta )

# Now calculate effect size for each group???????
model_NP <- rma(yi, vi, data = subset(meta, meta$climate_zone=="NP"))
summary.rma(model_NP)
#estimate      se     zval    pval    ci.lb   ci.ub 
# -0.0902    0.1086  -0.8307  0.4062  -0.3029  0.1226 


model_NT <- rma(yi, vi, data = subset(meta, meta$climate_zone=="NT"))
summary.rma(model_NT)
#estimate      se    zval    pval   ci.lb   ci.ub 
# 0.2693     0.0968  2.7826  0.0054  0.0796  0.4590  ** 

model_NST <- rma(yi, vi, data = subset(meta, meta$climate_zone=="NST"))
summary.rma(model_NST)
# estimate      se    zval    pval   ci.lb   ci.ub 
#   0.4506  0.0420  10.7289  <.0001  0.3683  0.5329  *** 

model_T <- rma(yi, vi, data = subset(meta, meta$climate_zone=="T"))
summary.rma(model_T)
#estimate      se    zval    pval   ci.lb   ci.ub 
#0.3370      0.0830  4.0584  <.0001  0.1742  0.4998  *** 

model_SST <- rma(yi, vi, data = subset(meta, meta$climate_zone=="SST"))
summary.rma(model_SST)
#estimate      se    zval    pval    ci.lb   ci.ub 
#0.2807     0.1597  1.7578  0.0788  -0.0323  0.5937  .

model_ST <- rma(yi, vi, data = subset(meta, meta$climate_zone=="ST"))
summary.rma(model_ST)
#estimate      se    zval    pval    ci.lb   ci.ub 
# 0.3740     0.2950  1.2676  0.2049  -0.2043  0.9522  

# Load these model outputs
model_output <- read_csv("Output/Data/model_outputs.csv") %>%
  mutate(climate_zone = factor(climate_zone, levels = c("ST", "SST", "T", "NST", "NT", "NP")))

# Plot?

ggplot(model_output, aes(y= climate_zone, x = effect_est, xmin = ci_lb, xmax = ci_ub)) +
  geom_point() +
  geom_errorbar(height = .1) +
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +
  theme_classic()


# another way of trying to stats
meta_model <- rma(yi, vi, data = meta) 
summary.rma(meta_model)
# so if we ignore everything, the overall effect of animals on pp is 0.3092, and this is significant
forest.rma(meta_model)


meta_model2 <- rma(yi, vi, mods = ~ climate_zone, data = meta )
summary.rma(meta_model2)
forest(meta_model2)

meta_model3 <- rma(yi, vi, mods = ~ lat_dec, data = meta )
summary.rma(meta_model3)
forest(meta_model3, cex = 1, subset(study))

