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

# mine data
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
  mutate(prov_N_treat = sum(nh4_treatment, no3_treatment, total_nitrogen_treatment, na.rm = TRUE),
         prov_N_cont = sum(nh4_control, no3_control, total_nitrogen_control, na.rm = TRUE),
         added_N = prov_N_treat - prov_N_cont)

# Graphing --------
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
ggplot(meta_final, aes(pp_taxa, climate_zone, colour = pp_taxa)) +
  geom_jitter(size = 4) +
  theme_black() +
  theme(legend.position = "none") +
  labs(x = "Primary Producer", y = "Climate Zone") +
  scale_y_discrete(labels = climate_labs) +
  scale_x_discrete(labels = c("Calcareous algae", "Coral", "Endosymbiont", "Fleshy Algae", "Mangrove", "Phytoplankton", "Seagrass")) +
  theme(axis.text.x = element_text(angle = 30))
  
# take a look at the distribution of consumers across latitudes
ggplot(meta_final, aes(c_taxa, climate_zone, colour = c_taxa)) +
  geom_jitter(size = 4) +
  theme_black() +
  theme(legend.position = "none") +
  labs(x = "Consumer", y = "Climate Zone") +
  scale_y_discrete(labels = climate_labs) +
  scale_x_discrete(labels = c("Birds", "Fish", "Crustaceans", "Echinoderms", "Gastropods", "Misc mobile inverts", "Bivalves", "Sponges", "Misc sessile inverts")) +
  theme(axis.text.x = element_text(angle = 30))

# look at consumer groups
ggplot(meta_final, aes(c_activity, climate_zone, colour = c_activity)) +
  geom_jitter(size = 4) +
  theme_black() +
  theme(legend.position = "none") +
  labs(x = "Consumer Activity Level", y = "Climate Zone") +
  scale_y_discrete(labels = climate_labs)

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

# Looking at background nutrients and temperature and provisioning ------

ggplot(meta_final, aes(abs_lat, background_nut)) +
  geom_point() +
  geom_smooth(method = lm)

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

# Latitude mulit level model ----

# We are going to use the absolute value of latitude, because I have no reason to expect North temperate to vary from South temperate

model_lat_paper <- rma.mv(yi = yi,
              V = vi,
              slab = author,
              data = meta_final,
              method = "REML",
              mods = ~ abs_lat,
              random = ~ 1 | paper_num/indiv, 
              test = "t")

summary(model_lat_paper)


# Look at Test of Moderators to determine whether there are significant differences between the subgroups
# p-value is 0.6181 so NO, there is NO difference!
# The Model Results are printed within a meta-regression framework. This means that we cannot directly extract the estimates in order to obtain the pooled effect sizes within subgroups.

# Graph lat ------

# Use predict() to predict effect size for each abs lat
# https://wviechtb.github.io/metafor/reference/predict.rma.html

# For random-effects models that were fitted with the rma.mv function, the model may actually include multiple 𝜏2 values (i.e., when the random argument includes an ‘~ inner | outer’ term and struct="HCS", struct="DIAG", struct="HAR", or struct="UN"). In that case, the function will provide prediction intervals for each level of the inner factor (since the prediction intervals differ depending on the 𝜏2 value)

# don't have to worry about tau2.levels and gamma2.levels arguments because I'm not using newmods, I'm just predicting the effects at the lat values from my input data

lat_predict <- as.data.frame(predict(model_lat_paper,
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

#ggsave("Output/Figures/abs_lat_effect", device = "png",
#     height = 9, width = 16, dpi = 400)

# Chlorophyll a multilevel models -----
# Ok there are multiple ways to look at chl a data

# Untransformed chl a estimates as a linear model
model_chla_paper_est <- rma.mv(yi = yi,
                               V = vi,
                               slab = author,
                               data = meta_final,
                               method = "REML",
                               mods = ~ chla_est,
                               random = ~ 1 | paper_num/indiv,
                               test = "t")

summary(model_chla_paper_est)
# p value is 0.5450

# Transformed log chl a data
# Do this when you know that a predictor wouldn't increase the response forever (non linear)
model_chla_paper_log <- rma.mv(yi = yi,
                               V = vi,
                               slab = author,
                               data = meta_final,
                               method = "REML",
                               mods = ~ chla_log,
                               random = ~ 1 | paper_num/indiv,
                               test = "t")
summary(model_chla_paper_log)
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


AIC(model_chla_paper_numeric, model_chla_paper_bin, model_chla_paper_log, model_chla_paper_est)
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
chla_predict <- as.data.frame(predict(model_chla_paper_est,
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

# Graph chl a log------
chla_predict_log <- as.data.frame(predict(model_chla_paper_log,
                                      addx=TRUE)) %>%
  tibble::rownames_to_column("author") %>%
  rename(chla_log = X.chla_log) 

# set seed because otherwise I'll get different colours each time
set.seed(1235)

# plot predicted values as line, jitter observed effect sizes
ggplot() +
  geom_line(data = chla_predict_log, aes(chla_log, pred)) +
  geom_jitter(data = meta_final, aes(chla_log, yi, colour = author)) +
  geom_ribbon(data = chla_predict_log, aes(chla_log, pred, ymin = ci.lb, ymax = ci.ub), alpha = .1) +
  theme_classic() +
  labs(x = "Log(chl a estimate)", y = "Predicted Effect Size") +
  theme(legend.position="none") +
  scale_color_manual(values = sample(col_vector, n)) +
  geom_hline(yintercept=0, color="black", linetype="dashed")

#ggsave("Output/Figures/chla_log_effect", device = "png",
 #      height = 9, width = 16, dpi = 400)

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

# Background nutrients multilevel model-----
model_nut_lat_paper <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_final,
                          method = "REML",
                          mods = ~ background_nut*abs_lat,
                          random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                          test = "t")

summary(model_back_paper)
# 0.3141

# Graph background nutrients ------
nut_predict <- as.data.frame(predict(model_back_paper,
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

#ggsave("Output/Figures/background_nut_effect", device = "png",
#     height = 9, width = 16, dpi = 400)

# Effect of temperature multilevel model-----
model_temp_paper <- rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = meta_final,
                           method = "REML",
                           mods = ~ temp * background_nut * abs_lat,
                           random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                           test = "t")

summary(model_temp_paper)
# Oh? temperature has an effect?

# check for colinearity
vif(model_temp_paper) # run this on the model without an interaction
# value of 1 means they're not correlated with any other variables
# higher value = greater correlation
# values of more than 4 or 5 are regarded as moderate - high
# values of 10+ are very high

# OK so there's some colinearity apparently 4.1667 is moderate

# Natural-ness of experiment multilevel model-----
model_nat_paper <- rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = meta_final,
                           method = "REML",
                           mods = ~ natural,
                           random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                           test = "t")

summary(model_nat_paper)
# Look at Test of Moderators to detemine whether there are significant differences between the subgroups
# p-value is 0.3003 so no, there is no difference!
# The Model Results are printed within a meta-regression framework. This means that we cannot directly extract the estimates in order to obtain the pooled effect sizes within subgroups.

# How many papers in each bin?
nat_n <- meta_final %>%
  group_by(natural) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:5,
         nat_bin = factor(c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm"),
                          levels = c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")))

# Pull coefficient estimates so I can graph them
intercept_nat <- coef(summary(model_nat_paper))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes

# edit bin labels
nat_effect_sum <- coef(summary(model_nat_paper)) %>%
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

nat_n_labs <- c("n = 52", "n = 14", "n = 16", "n = 9", "n = 43")

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

ggsave("Output/Figures/effects_nat", device = "png",
     height = 9, width = 8, dpi = 400)

# Consumer activity level multilevel model -----
model_act_paper <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_final,
                          method = "REML",
                          mods = ~ c_activity,
                          random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                          test = "t")

summary(model_act_paper)
# Look at Test of Moderators to detemine whether there are significant differences between the subgroups
# p-value is 0.4606 so no, there is no difference!
# The Model Results are printed within a meta-regression framework. This means that we cannot directly extract the estimates in order to obtain the pooled effect sizes within subgroups.

# How many papers in each bin?
act_n <- meta_final %>%
  group_by(c_activity) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:3)

# Pull coefficient estimates so I can graph them
intercept_act <- coef(summary(model_act_paper))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes

act_effect_sum <- coef(summary(model_act_paper)) %>%
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

act_n_labs <- c("n = 44", "n = 21", "n = 69")

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

ggsave("Output/Figures/effects_act", device = "png",
     height = 9, width = 8, dpi = 400)


# I feel like at some point I will have to address the differences between the tropical taxons and the temperate taxons


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

# Separate climate zones -------
tropical <- meta_final %>%
  filter(climate2 == "Tropics")

subtropics <- meta_final %>%
  filter(climate2 == "Subtropics")

temperate <- meta_final %>%
  filter(climate2 == "Temperate")

# Does chl a matter in the tropics?
model_chla_tropical <- rma.mv(yi = yi,
                               V = vi,
                               slab = author,
                               data = tropical,
                               method = "REML",
                               mods = ~ chla_est,
                               random = ~ 1 | paper_num/indiv,
                               test = "t")

summary(model_chla_tropical) #No

# Does chl a matter in the subtropics?
model_chla_subtropics <- rma.mv(yi = yi,
                              V = vi,
                              slab = author,
                              data = subtropics,
                              method = "REML",
                              mods = ~ chla_est,
                              random = ~ 1 | paper_num/indiv,
                              test = "t")

summary(model_chla_subtropics) #Not really

# Does chl a matter in the temperates?
model_chla_temperate <- rma.mv(yi = yi,
                              V = vi,
                              slab = author,
                              data = temperate,
                              method = "REML",
                              mods = ~ chla_est,
                              random = ~ 1 | paper_num/indiv,
                              test = "t")

summary(model_chla_temperate) #No

# Don't bother with polar bc only one study

# Independence and Bias -------

funnel.rma(model_lat_paper)
funnel.rma(model_chla_paper_est)
funnel.rma(model_nat_paper)
# A bit of a right skew

fsn(yi = yi, vi = vi, data = meta_final, type="Rosenthal") #file drawer analysis (missing studies have an effect size = 0)
# Fail-safe N: 44,981
# That seems like quite a large number?????

# Check for influential studies? ------

# make a df where each row is a paper, and num_obs is the number of observations that paper contributed
max_obs <- meta_final %>%
  group_by(paper_num) %>%
  dplyr::summarise(num_obs = as.numeric(n()))

hist(max_obs$num_obs)

total_num_obs <- sum(max_obs$num_obs)

# there's 173 data points
# so 10% of that would be 17.3 observations per paper
# so that would mean that paper # 17 and #23 and #29 have too many. try running things without those two

# Ok so let's remove those papers one at a time and rerun the stats???

meta_drop <- meta_final %>%
  filter(paper_num != 29) 

# removing 31 really changes the lat result
  # contributed the highest 2 effect sizes
# and removing 23 really changes chl a result
  # 23 was the highest chl a estimate
# see what happens when I remove both

# Rerun lat ----
model_lat_drop <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_drop,
                          method = "REML",
                          mods = ~ abs_lat,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_lat_drop)

# Removed # 17, p-value is 0.7132 versus 0.6181
# Removed # 23, p-value is 0.6927 versus 0.6181
# Removed # 31, p-value is 0.3172 versus 0.6181

# Nothing really changes things too much

# lets graph that
lat_predict_drop <- as.data.frame(predict(model_lat_drop,
                                     addx=TRUE)) %>%
  tibble::rownames_to_column("author") %>%
  rename(abs_lat = X.abs_lat) 

# set seed because otherwise I'll get different colours each time
set.seed(1235)

# plot predicted values as line, jitter observed effect sizes
ggplot() +
  geom_line(data = lat_predict_drop, aes(abs_lat, pred)) +
  geom_jitter(data = meta_drop, aes(abs_lat, yi, colour = author)) +
  geom_ribbon(data = lat_predict_drop, aes(abs_lat, pred, ymin = ci.lb, ymax = ci.ub), alpha = .1) +
  theme_classic() +
  labs(x = "Absolute latitude", y = "Predicted Effect Size") +
  theme(legend.position="none") +
  scale_color_manual(values = sample(col_vector, n)) +
  geom_hline(yintercept=0, color="black", linetype="dashed")


# Rerun chl a -----
meta_drop <- meta_final %>%
  filter(paper_num != 23)

meta_23 <- meta_final %>%
  filter(paper_num == 23)

meta_dropped <- meta_final %>%
  filter(paper_num == 31) %>%
  rbind(meta_23) %>%
  select(author, lat_dec, season,year_of_study, abs_lat, climate_zone, yi, vi, chla_est)

model_chla_drop <- rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = meta_drop,
                           method = "REML",
                           mods = ~ chla_est,
                           random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                           test = "t")

summary(model_chla_drop)
# Remove # 17 p-value is .6768 versus 0.5450 
# Remove # 23 p-value is 0.0806 versus 0.5450 # OK BIG CHANGE
# Remove # 29 p-value is 0.4960 versus 0.5450 


chla_predict_drop <- as.data.frame(predict(model_chla_drop,
                                      addx=TRUE)) %>%
  tibble::rownames_to_column("author") %>%
  rename(chla_est = X.chla_est) 

# set seed because otherwise I'll get different colours each time
set.seed(1235)

# plot predicted values as line, jitter observed effect sizes
ggplot() +
  geom_line(data = chla_predict_drop, aes(chla_est, pred)) +
  geom_jitter(data = meta_drop, aes(chla_est, yi, colour = author)) +
  geom_ribbon(data = chla_predict_drop, aes(chla_est, pred, ymin = ci.lb, ymax = ci.ub), alpha = .1) +
  theme_classic() +
  labs(x = "Chl a estimate", y = "Predicted Effect Size") +
  theme(legend.position="none") +
  scale_color_manual(values = sample(col_vector, n)) +
  geom_hline(yintercept=0, color="black", linetype="dashed")

# OK so I think I'll have to drop # 23 and # 31 
# And that means both lat and chl a are not significant predictors

# Rerun interaction model ------
meta_drop <- meta_final %>%
  filter(paper_num != 23)

model_multi_drop <- rma.mv(yi = yi,
                      V = vi,
                      slab = author,
                      data = meta_drop,
                      method = "REML",
                      mods = ~ chla_est * abs_lat,
                      random = ~ 1 | paper_num/indiv,
                      test = "t")

summary(model_multi_drop)

# check for colinearity
vif(model_multi_drop) # run this on the model without an interaction

# Remove # 17 p-value is 0.0041 versus p-val = 0.0194 
# Remove # 23 p-value is 0.0970 versus p-val = 0.0194
  # but the p val for the chl a, lat, and interaction are still basically the same
# Remove # 29 p-value is 0.0241 versus p-val = 0.0194 

# OK so it looks like I might not need to drop #23 ??????

# Rerun background nutrients multilevel model-----

meta_drop <- meta_final %>%
  filter(paper_num != 29)

model_back_paper <- rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = meta_drop,
                           method = "REML",
                           mods = ~ background_nut,
                           random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                           test = "t")

summary(model_back_paper)
# 0.3141
# Remove # 17 p-value is 0.2242 versus p-val =  0.3141
# Remove # 23 p-value is 0.1181 versus p-val =  0.3141
# Remove # 29 p-value is 0.7017 versus p-val =  0.3141

# Rerun natural-ness of experiment -----
re_model_nat_paper <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_f,
                          method = "REML",
                          mods = ~ natural,
                          random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                          test = "t")

summary(re_model_nat_paper)

# Remove # 17 p value is 0.2500 vs 0.3003 
# Remove # 23 p value is 0.3453 vs 0.3003 
# Remove # 23 p value is 0.3003 vs 0.1630 


# How many papers in each bin?
re_nat_n <- meta_f %>%
  group_by(natural) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:5,
         nat_bin = factor(c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm"),
                          levels = c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")))

# Pull coefficient estimates so I can graph them
re_intercept_nat <- coef(summary(re_model_nat_paper))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes

# edit bin labels
re_nat_effect_sum <- coef(summary(re_model_nat_paper)) %>%
  mutate(nat_bin = factor(c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm"),
                          levels = c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")),
         effect = ifelse(nat_bin != "Field: unstocked", estimate + re_intercept_nat,
                         re_intercept_nat),
         ci.lb = ifelse(nat_bin != "Field: unstocked", ci.lb + re_intercept_nat,
                        ci.lb),
         ci.ub = ifelse(nat_bin != "Field: unstocked", ci.ub + re_intercept_nat,
                        ci.ub)) %>% # add intercept to other coefficients
  left_join(re_nat_n, by = "nat_bin")

# GRAPH THIS?
nat_labs <- c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")

re_nat_n_labs <- c("n = 72", "n = 14", "n = 16", "n = 10", "n = 44")

ggplot(re_nat_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = nat_bin)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size (ROM)", y = " ") +
  scale_y_continuous(breaks = 1:length(nat_labs),
                     labels = nat_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(re_nat_n_labs),
                                         labels = re_nat_n_labs)) +
  scale_colour_manual(values = rev(pal_star)) +
  theme_black()+ 
  theme(legend.position = "none")

# Removed # 17 basically the same, field unstocked now barely crosses 0
# Removed # 23 basically the same
# Removed # 31 basically the same

# Rerun consumer activity level -----
re_model_act_paper <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_f,
                          method = "REML",
                          mods = ~ c_activity,
                          random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                          test = "t")

summary(re_model_act_paper)
# Removed # 17 p-value is 0.5830 vs 0.4606 
# Removed # 23 p-value is 0.4606 vs 0.4606 
# Removed # 31 p-value is 0.1156 vs 0.4606 


# How many papers in each bin?
re_act_n <- meta_f %>%
  group_by(c_activity) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:3)

# Pull coefficient estimates so I can graph them
re_intercept_act <- coef(summary(re_model_act_paper))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes

re_act_effect_sum <- coef(summary(re_model_act_paper)) %>%
  mutate(c_activity = factor(c("Vertebrate", "Mobile Invert", "Sessile Invert"),
                             levels = c("Vertebrate", "Mobile Invert", "Sessile Invert")),
         effect = ifelse(c_activity != "Vertebrate", estimate + re_intercept_act,
                         re_intercept_act),
         ci.lb = ifelse(c_activity != "Vertebrate", ci.lb + re_intercept_act,
                        ci.lb),
         ci.ub = ifelse(c_activity != "Vertebrate", ci.ub + re_intercept_act,
                        ci.ub)) %>% # add intercept to other coefficients
  left_join(re_act_n, by = "c_activity")

# GRAPH THIS?
act_labs <- c("Vertebrate", "Mobile Invert", "Sessile Invert")

re_act_n_labs <- c("n = 46", "n = 21", "n = 89")

ggplot(re_act_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = c_activity)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size (ROM)", y = " ") +
  scale_y_continuous(breaks = 1:length(act_labs),
                     labels = act_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(re_act_n_labs),
                                         labels = re_act_n_labs)) +
  scale_colour_manual(values = rev(pal_star)) +
  theme_black()+ 
  theme(legend.position = "none")

# Same when I remove # 17 and 23 and 31

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



# CODE GRAVEYARD ---------

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
intercept_lat <- coef(summary(model_lat_paper))[1,1] # extract the intercept

# How many papers in each bin?
lat_n <- meta_final %>%
  group_by(climate_zone) %>%
  summarise(n = n()) %>%
  rename(climate = climate_zone) %>%
  mutate(dummy = 1:6)

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes
lat_effect_sum <- coef(summary(model_lat_paper)) %>%
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

