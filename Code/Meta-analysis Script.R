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
library(patchwork)

source(file = "Code/theme_black_code.R")
theme_set(theme_classic())

# Set preferences so I can see all my columns
rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

# mine data from papers using this site
# https://automeris.io/WebPlotDigitizer/


# Load data
# And do a million manipulations......
meta_full <- read_csv("Output/Data/march_meta2023.csv") %>%
  mutate(sd_treatment = replace_na(sd_treatment, 999), # assigning a dummy to replace
         sd_control = replace_na(sd_control, 999)) %>% # assigning a dummy to replace
  mutate(sd_treatment = ifelse(sd_treatment == 999, se_treatment*sqrt(treatment_n), sd_treatment), # replacing that dummy with the estimated SD
         sd_control = ifelse(sd_control == 999, se_control*sqrt(control_n), sd_control),# replacing that dummy with the estimated SE
         natural = as.factor(natural),
         pp_taxa = as.factor(pp_taxa),
         abs_lat = abs(lat_dec), # creating an absolute value lat
         climate = factor(case_when(abs_lat <= 23.45 ~ "Tropics",
                                     abs_lat <= 35 ~ "Subtropics",
                                     abs_lat <= 60 ~ "Temperate",
                                     abs_lat > 60 ~ "Polar"),
         levels = c("Tropics", "Subtropics", "Temperate", "Polar")), # New column where I don't distinguish between North subtropics/temperates and South subtropics/temperates
         c_group = case_when(c_taxa == "Bi" ~ "Bird",
                             c_taxa == "F" ~ "Fish",
                             c_taxa == "E" ~ "Echinoderm",
                             c_taxa == "G" ~ "Mollusc", 
                             c_taxa == "B" ~ "Mollusc",
                             c_taxa == "SP" ~ "Sponge",
                             c_taxa == "SI" ~ "Cnidaria",
                             c_taxa == "MI" ~ "Cnidaria"),
         c_group = factor(c_group, levels = c("Bird", "Fish", "Echinoderm", "Mollusc", "Cnidaria", "Sponge")),
         nat_bin = case_when(natural == 1 ~ "Field: unstocked",
                             natural == 2 ~ "Field: PP stocked",
                             natural == 3 ~ "Field: animals stocked",
                             natural == 4 ~ "Field: both stocked",
                             natural == 5 ~ "Mesocosm"),
         nat_bin = factor(nat_bin, levels = c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")))

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

meta_calc <- meta_no_neg %>%
  escalc(measure = "ROM", data = meta_no_neg, # ROM for log transformed ratio of means
       m1i = mean_treatment, sd1i = sd_treatment, # log(treatment/control)
       m2i = mean_control, sd2i = sd_control,
       n1i = treatment_n, n2i = control_n) %>%
  rowwise() %>%
  mutate(sum_all_N_treat = sum(nh4_treatment, no3_treatment, total_nitrogen_treatment, na.rm = TRUE),
         sum_all_N_cont = sum(nh4_control, no3_control, total_nitrogen_control, na.rm = TRUE),
         added_N = na_if((sum_all_N_treat - sum_all_N_cont), 0),
         added_nh4 = na_if((sum(nh4_treatment, na.rm = TRUE) - sum(nh4_control, na.rm = TRUE)), 0)) %>%
  ungroup() %>%
  mutate(pp_group = as.factor(case_when(
          pp_taxa == "CA" | pp_taxa == "FA" ~ 'Algae',
          pp_taxa == "En" | pp_taxa == "Co" ~ 'Coral',
          pp_taxa == "P" ~ 'Phytoplankton',
          pp_taxa == "S" ~ "Plants",
          pp_phyla == "Tracheophyta" ~ "Plants",
          pp_taxa == "O" & pp_phyla == "Cnidaria" ~ "Coral")),
         pp_groupz = as.factor(case_when(
           pp_taxa == "CA"  ~ 'Crustose algae',
           pp_taxa == "FA" ~ 'Fleshy algae',
           pp_taxa == "En" ~ "Endo",
           pp_taxa == "Co" ~ 'Coral',
           pp_taxa == "P" ~ 'Phytoplankton',
           pp_taxa == "S" ~ "Seagrass",
           pp_taxa == "O" & pp_phyla == "Tracheophyta" ~ "Trees",
           pp_taxa == "O" & pp_phyla == "Cnidaria" ~ "Coral")),
         lat_stand = base::scale(abs_lat),
         background_nut_stand = scale(background_nut),
         added_N_stand = scale(added_N),
         temp_stand = scale(temp))

# Want all models to use the same data!
# this df only makes sure to standardize the df between the three correlated variables, and is for use in models without added_N
meta_data <- meta_calc %>%
  drop_na(background_nut_stand) %>%
  drop_na(temp_stand) %>%
  mutate(indiv = as.factor(indiv),
         paper_num = as.factor(paper_num))

# Multilevel models ------
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
# If this p-value > 0.05 then we can reject null hypothesis H0 : β1 = β2 = 0 based on the omnibus test (this means our moderators are better than a null model)

# The Model Results are printed within a meta-regression framework. This means that we cannot directly extract the estimates in order to obtain the pooled effect sizes within subgroups.

# Moderators: -----
# Main question: is effect size altered by latitude?
# Main mods: 
  # temperature
  # latitude
  # background nutrients

# Co-mods:
  # provisioned nutrients
  # openness of experiment
  # PP group (algae, coral, phytoplankton or plant)

# Why did I stop looking at consumer activity level?

# When I build a separate model for each moderator the coefficients and p-values don't change
# Keep them in one model!

# A note on the added_N 
# I built models that included added_N and it didn't improve the model or have any effect and cuts like 100 rows from my data so I've decided to drop it
# The model with added_N and a smaller sample size comes out looking EXACTLY like the model without added_N and a larger sample size. 


# One model to rule them all ------
# If goal is to explain patterns, and sample size is big enough to include all the predictors, then just throw them all in the same model. 
# Standardize the numeric predictors, then compare standardized effect sizes
# AIC is good for when you can’t put everything in the same model, I might not be able to put temperature and lat in the same model

model_all <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_data,
                          method = "REML",
                          mods = ~ lat_stand + background_nut_stand  + temp_stand + natural + pp_group + lat_stand:background_nut_stand +  temp_stand:background_nut_stand ,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_all)

metafor::vif(model_all) # Lat and temp are correlated
# intercept tells me what the value is at the natural 1 pp group algae and mean of all other variables
# looking at the full data (took out all added_N terms) the coefficients are basically hte same
# so we can look at the reduced data and it's saying something very similar to the full data

# In general, all three ways I’ve tried to do these stats give a similar answer. Identity of the pp matters (algae has a stronger response than others), and background nutrients is weak/moderate positive effect. This makes both # Hannah and I feel better. All methods point to context dependence!

# Correlation between variables
cor(cbind(meta_data$lat_stand, meta_data$background_nut_stand, meta_data$temp_stand), use = "complete.obs")
# so temp, lat, and background_nutrients are correlated (|value| < 0.5)
# So maybe we should go back to AIC because we shouldn't have these three in the same model
# BC they're correlated, there might be a generally negative relationship between lat/temp/back_nut and effect but within that slope the model things back_nut has a positive effect


# Let's make sure the three level model is better than a two level model before I go further
# fit a model in which the level 3 variance, representing the between-study heterogeneity, is set to zero
# This is equal to fitting a simple random-effects model in which we assume that all effect sizes are independent (which we know they are not).

model_all_l3_removed <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = meta_data,
                    method = "REML",
                    mods = ~ lat_stand + background_nut_stand  + temp_stand + natural + pp_group + lat_stand:background_nut_stand +  temp_stand:background_nut_stand,
                    random = ~ 1 | paper_num/indiv, 
                    test = "t",
                    sigma2 =  c(0, NA)) # set level 3 variance = 0

summary(model_all_l3_removed) # We can see we sucessfully set sigma 2 (level 3 variance) = 0

# But is it better?
anova(model_all, model_all_l3_removed)
# They're basically the same
# BUT the website says even if the three level model isn't better, if we KNOW the studies are not independent it still makes sense to model the random effect of study to be statistically correct
# So I think I'm justified in keeping the third level
# And the effect sizes and p-values don't change in a meaningful way, if at all.


#### Three separate models ------
# 1) everything in full model - temp and lat
# 2) everything in full model - temp and back_nut
# 3) everything in full model - back_nut and lat

# Then compete them using AIC!

# rational for AIC is that we can't put all three in the same model bc they mess up each other
# AIC will tell us which of the three is best

#Subset each plant type to rerun the model_all with each plant type? Or just to do algae by itself and then the rest by themselves. 
# -might have to run each mod by itself bc low sample size

# Except I'm going to compete the "full" model with lat vs "full" with temp vs "full" with background_nutrients


# full model with just latitude
model_full_lat <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = meta_data,
                    method = "REML",
                    mods = ~ lat_stand  + natural + pp_group,
                    random = ~ 1 | paper_num/indiv, 
                    test = "t")

summary(model_full_lat)

# Full model with background_nutrients
model_full_back <- rma.mv(yi = yi,
                                V = vi,
                                slab = author,
                                data = meta_data,
                                method = "REML",
                                mods = ~ background_nut_stand  + natural + pp_group,
                                random = ~ 1 | paper_num/indiv, 
                                test = "t")

summary(model_full_back)


# Full model with temperature
model_full_temp <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_data,
                          method = "REML",
                          mods = ~ temp_stand + natural + pp_group,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_full_temp)

# Null model
model_null <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_data,
                          method = "REML",
                          mods = ~ 1,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_null)

# Model with just pp and nat
model_cat <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_data,
                          method = "REML",
                          mods = ~ pp_group + natural,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_cat)

# Now see which model has the best support
AIC_tab <- AIC(model_full_lat, model_full_back, model_full_temp, model_null, model_cat) %>%
  rownames_to_column() %>%
  mutate(best = min(AIC),
         delta = AIC - best)

# So the model with background nutrients is the best model according to AIC but there's still moderate support for all three of lat, temp, and background nutrients

# Present the coefficient plot of the best model, but include the other two in the supp
# And we still get an estimate for the effect size for everything
# Since none of the models are better than the other, I can just chose the one I'm most interested in and throw the other two in the supp
# But.... I'm interested in both lat and background nut

# AIC tells you which collection of variables best explain the data
# And then looking at the p-values of those variables tells you whether ___
# Hannah advocates for just looking at the effect size (coefficient)


# Figures --------
# Maybe some help here?
# https://www.metafor-project.org/doku.php/tips:meta_regression_with_log_rr

# Define a palette
pal_star <- pnw_palette("Starfish", 5)
pal_sm <- pnw_palette("Starfish", 3)

# Calculations for plot effect sizes -------
# Background nutrients
# Pull coefficient estimates so I can graph them
back_nut_effect_sum <- coef(summary(model_full_back)) %>%
  slice(-(3:9)) %>%
  slice(-1)%>%
  rownames_to_column()%>%
  mutate(dummy = 1, effect = estimate)

# Experimental design
nat_n <- meta_data %>%
  group_by(nat_bin) %>%
  summarise(n = n()) %>% # How many papers in each bin?
  mutate(dummy = 1:5)

# Just the intercept to adjust other estimates
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


# Primary producer group
pp_n <- meta_data %>%
  group_by(pp_group) %>%
  summarise(n = n()) %>% # How many papers in each bin?
  mutate(dummy = 1:4)

# Just the intercept to adjust other estimates
intercept_pp <- coef(summary(model_full_back))[1,1] # extract the intercept

# Pull rest of the coefficients, add a column for pp, and adjust the estimates to reflect the estimated effect sizes
pp_effect_sum <- coef(summary(model_full_back)) %>%
  slice(-(2:6)) %>%
  mutate(pp_group = factor(c("Algae", "Coral", "Phytoplankton", "Plants")), 
         effect = ifelse(pp_group != "Algae", estimate + intercept_pp,
                         intercept_pp), # add intercept to the estimated coefficient
         ci.lb_adjust = effect - (1.96*se), #adjust the CI based on the new estimate
         ci.up_adjust = effect + (1.96*se)) %>% 
  left_join(pp_n, by = "pp_group")

# Labels
cont_labs <- c("Background nutrients", "Absolute latitude", "Temperature")
names(cont_labs) <- c("background_nut_stand", "lat_stand", "temp_stand")
back_n_labs <- c("n = 172")

nat_labs <- c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")
nat_n_labs <- c("n = 74", "n = 14", "n = 21", "n = 9", "n = 54")

pp_labs <- c("Algae", "Coral", "Phytoplankton", "Plants")
pp_n_labs <- c("n = 25", "n = 47", "n = 53", "n = 47")

# plot each moderator using the function I created (plot_themes)
back_nut_plot <- effect_plot(back_nut_effect_sum, ci.lb, ci.ub, rowname, pal_sm[3], cont_labs, back_n_labs) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())

nat_plot <- effect_plot(nat_effect_sum, ci.lb_adjust,ci.up_adjust, nat_bin, rev(pal_star), nat_labs, nat_n_labs) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) 

pp_plot <- effect_plot(pp_effect_sum, ci.lb_adjust, ci.up_adjust, pp_group, rev(pal_star), pp_labs, pp_n_labs) +
  xlab("Effect Size")

# Effect size plot ------
back_nut_plot + nat_plot + pp_plot + 
  plot_layout(ncol = 1, heights = unit(c(1, 4, 5), c("cm", "null"))) + plot_annotation(tag_levels =  "A")

#ggsave("Output/Pub_figs/full_back_nut_mod.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# Graph the other models!
# Lat effect size
lat_effect_sum <- coef(summary(model_full_lat)) %>%
  slice(-(3:9)) %>%
  slice(-1) %>%
  rownames_to_column() %>%
  mutate(dummy = 1, effect = estimate)


# Temperature effect size
temp_effect_sum <- coef(summary(model_full_temp)) %>%
  slice(-(3:9)) %>%
  slice(-1) %>%
  rownames_to_column() %>%
  mutate(dummy = 1, effect = estimate)

# Consumers vs pp vs lat
pal <- pnw_palette("Bay", 6)
pal2 <- pnw_palette("Bay", 4)

meta_data %>% 
  mutate(pp_taxa = factor(pp_taxa, levels = rev(c("Co", "En", "P", "CA", "FA", "S", "O")))) %>%
  ggplot(aes(c_group, pp_taxa, colour = climate, fill = climate, pch = nat_bin)) +
  geom_jitter(size = 6, alpha = 0.3, width = 0.25, height = 0.25) +
  theme_white() +
  labs(x = "Consumer taxa", y = "Primary producer taxa", colour = "Climate",
       pch = "Experiment", fill = "Climate") +
  theme(axis.text.x = element_text(angle = 30)) +
  scale_y_discrete(labels = c("Mangrove", "Seagrass", "Fleshy Algae", "Calcareous algae", "Phytoplankton", "Endosymbiont", "Coral")) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)) +
  scale_colour_manual(values = rev(pal2)) +
  scale_fill_manual(values = rev(pal2)) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25))

#ggsave("Output/Pub_figs/taxa_region.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# Tables ------

# Let's make a table of the different pp responses! -----
# OK SO what I think I want is twofold
# I want to present a table of the metrics I got means of to estimate productivity
# That will be something like growth (change/time), biomass, percent cover, density, etc
# AND THEN I probably want to report what sorts of things people were measuring the growth rates of!
# I went through and did a bunch of post - pre calculations for studies where they gave me measurements over timepoints, this counts as growth

metrics <- meta_data %>%
  mutate(response = case_when(response == "mass" ~ "biomass",
                              TRUE ~ as.character(response))) %>%
  count(response)

# Look at individual metrics to sea what was measured for each
measured <- meta_data %>%
  filter(response == "width") %>%
  select(response, measured, growth_rate, pp_taxa, comments)

# Which metrics were used to measure growth rate?
growth_metrics <- meta_data %>%
  filter(growth_rate == "Y") %>%
  count()

# What metrics were used when growth wasn't measured
other_metrics <- meta_data %>%
  filter(growth_rate == "N") %>%
  count(response)

# Which species were used? -----
pp_used <- meta_data %>%
  count(pp_group, pp_family, pp_genus, pp_species)

c_used <- meta_data %>%
  count(c_group, c_family, c_genus, c_species)

# Independence and Bias -------

funnel.rma(model_all)
funnel.rma(model_full_lat)
funnel.rma(model_full_back)
funnel.rma(model_full_temp)
funnel.rma(model_null)

# A bit of a right skew

fsn(yi = yi, vi = vi, data = meta_data, type="Rosenthal") #file drawer analysis (missing studies have an effect size = 0)
# Fail-safe N: 42,894
# That seems like quite a large number?????

# Check for influential studies? ------
# maybe I can use metafor::influence() to help
# ?metafor::influence.rma.mv()
# These tests take FOREVER don't run!

# cooks.distance(model_full_back, progbar=FALSE, reestimate=TRUE, parallel="no", ncpus=1, cl=NULL)

# dfbetas(model_full_back, progbar=FALSE, reestimate=TRUE, parallel="no", ncpus=1, cl=NULL, ...)

# hatvalues(model_full_back, type="diagonal")

# Do this manually instead
# make a df where each row is a paper, and num_obs is the number of observations that paper contributed
max_obs <- meta_data %>%
  group_by(paper_num) %>%  
  summarize(n_obs = n_distinct(obs_num))

hist(max_obs$n_obs)

sum(max_obs$n_obs) # 170
# so 10% of that would be 17 observations per paper
# so that would mean that paper # 17 and #23 and #29 have just barely too many

# Ok so let's remove those papers one at a time and rerun the stats! let's just use the one model bc fuck that

# # 17 and #23 and #29 are the three I need to check
meta_data_drop <- meta_data %>%
  filter(paper_num != 17) 

# Look at the full model
model_all_drop <- rma.mv(yi = yi,
                    V = vi,
                    slab = author,
                    data = meta_data_drop,
                    method = "REML",
                    mods = ~ lat_stand + background_nut_stand  + temp_stand + natural + pp_group + lat_stand:background_nut_stand +  temp_stand:background_nut_stand,
                    random = ~ 1 | paper_num/indiv, 
                    test = "t")
summary(model_all_drop)

# Summary stats from full model with all papers:
#Test of Moderators (coefficients 2:13):
#  F(df1 = 12, df2 = 159) = 1.9483, p-val = 0.0324
#
#Model Results:
#  
#                                estimate      se     tval   df    pval    ci.lb   ci.ub 
#intrcpt                            0.4731  0.1885   2.5102  159  0.0131   0.1009  0.8453  * 
#  lat_stand                         -0.1945  0.1412  -1.3778  159  0.1702  -0.4733  0.0843    
#background_nut_stand               0.1233  0.1458   0.8461  159  0.3988  -0.1645  0.4112    
#temp_stand                        -0.1165  0.1241  -0.9390  159  0.3492  -0.3615  0.1285    
#natural2                          -0.1699  0.2428  -0.6999  159  0.4850  -0.6494  0.3096    
#natural3                           0.3753  0.2025   1.8531  159  0.0657  -0.0247  0.7753  . 
#natural4                           0.5523  0.3119   1.7708  159  0.0785  -0.0637  1.1683  . 
#natural5                           0.1550  0.2331   0.6652  159  0.5069  -0.3053  0.6154    
#pp_groupCoral                     -0.1655  0.2170  -0.7624  159  0.4469  -0.5941  0.2631    
#pp_groupPhytoplankton             -0.4867  0.2489  -1.9553  159  0.0523  -0.9783  0.0049  . 
#pp_groupPlants                    -0.4024  0.2103  -1.9131  159  0.0575  -0.8178  0.0130  . 
#lat_stand:background_nut_stand    -0.1583  0.1595  -0.9923  159  0.3226  -0.4733  0.1567    
#background_nut_stand:temp_stand   -0.2366  0.1920  -1.2325  159  0.2196  -0.6157  0.1425  

# Removing # 17 makes latitude negative and significant in the full model
# Removing # 23 and 29 doesn't do anything

# Rerun with lat model
model_full_lat_drop <- rma.mv(yi = yi,
                               V = vi,
                               slab = author,
                               data = meta_data_drop,
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
                          data = meta_data_drop,
                          method = "REML",
                          mods = ~ background_nut_stand  + natural + pp_group,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_full_back_drop)

# Original model output: 
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

# Rerun with with temperature model
model_full_temp_drop <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_data_drop,
                          method = "REML",
                          mods = ~ temp_stand + natural + pp_group,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(model_full_temp_drop)

# Original model output: 
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

# Summary statistics! -----
# what percentage of studies fell into each category?
# Experimental design
nat_n2 <- nat_n %>%
  mutate(percent = n/sum(n),
         delta = 0.2 - percent)

# PP group
pp_n2 <- pp_n %>%
  mutate(percent = n/sum(n))

# Abs lat
lat_n <- meta_data %>%
  distinct(paper_num, .keep_all = TRUE) %>%
  count(climate)%>%
  mutate(percent = n/sum(n))

#plot distribution of latitude
meta_data %>%
  distinct(paper_num, .keep_all = TRUE) %>%
  mutate(dummy = 1) %>%
  ggplot(aes(abs_lat, dummy)) +
  geom_point()