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

# Packages for map
library(sf)
library(rnaturalearth) # another simple world map?


source(file = "Code/Plot_themes.R")
theme_set(theme_classic())

# Set preferences so I can see all my columns
rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

# mine data from papers using this site
# https://automeris.io/WebPlotDigitizer/

# Load data
# And do a million manipulations......
meta_full <- read_csv("Output/Data/april_meta2023.csv") %>%
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
                             c_taxa == "SI" ~ "Cnidaria", # this is the 1 obelia
                             c_taxa == "MI" ~ "Cnidaria"), # these are the jellies
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
  escalc(measure = "ROM", data = meta_no_neg, 
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
          pp_taxa == "O" & pp_phyla == "Cnidaria" ~ "Coral")), # the one anemone
         pp_groupz = as.factor(case_when(
           pp_taxa == "CA"  ~ 'Calcareous algae',
           pp_taxa == "FA" ~ 'Fleshy algae',
           pp_taxa == "En" ~ "Endosymbiont",
           pp_taxa == "Co" ~ 'Coral',
           pp_taxa == "P" ~ 'Phytoplankton',
           pp_taxa == "S" ~ "Seagrass",
           pp_taxa == "O" & pp_phyla == "Tracheophyta" ~ "Mangrove+cordgrass",
           pp_taxa == "O" & pp_phyla == "Cnidaria" ~ "Coral")),
         pp_groupz = factor(pp_groupz, levels = rev(c("Coral", "Endosymbiont", "Phytoplankton", "Calcareous algae", "Fleshy algae", "Seagrass", "Mangrove+cordgrass"))),
         lat_stand = base::scale(abs_lat),
         background_nut_stand = scale(background_nut_umolL),
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
AIC_tab <- AIC(model_full_lat, model_full_back, model_full_temp, model_null) %>%
  rownames_to_column() %>%
  mutate(best = min(AIC),
         delta = AIC - best,
         likelihood = exp( -0.5*delta),
         sum = sum(likelihood),
         AICw = likelihood/sum)

aic_table(model_full_lat, model_full_back, model_full_temp, model_null, model_cat)
# So the model with background nutrients is the best model according to AIC but there's still moderate support for all three of lat, temp, and background nutrients

# The SMD models have almost double the AIC value of the ROM models

# Present the coefficient plot of the best model, but include the other two in the supp
# And we still get an estimate for the effect size for everything
# Since none of the models are better than the other, I can just chose the one I'm most interested in and throw the other two in the supp
# But.... I'm interested in both lat and background nut

# AIC tells you which collection of variables best explain the data
# And then looking at the p-values of those variables tells you whether ___
# Hannah advocates for just looking at the effect size (coefficient)


# Figures --------

# fuck around and see what those new slattery papers did
meta_data %>%
  mutate(slat_new = ifelse(paper_num == 17 & obs_num > 15, "new_slat", "not")) %>%
ggplot(aes(background_nut_umolL, yi, colour = slat_new)) +
  geom_point()

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
         ci.up_adjust = effect + (1.96*se),
         duh = "yes") %>% 
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
         ci.up_adjust = effect + (1.96*se),
         duh = "yes") %>% 
  left_join(pp_n, by = "pp_group")

# Labels
cont_labs <- c("Background nutrients", "Absolute latitude", "Temperature")
names(cont_labs) <- c("background_nut_stand", "lat_stand", "temp_stand")
back_n_labs <- c("n = 180")

nat_labs <- c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")
nat_n_labs <- c("n = 84", "n = 12", "n = 21", "n = 9", "n = 54")

pp_labs <- c("Algae", "Coral", "Phytoplankton", "Plants")
pp_n_labs <- c("n = 27", "n = 55", "n = 53", "n = 45")

# plot each moderator using the function I created (plot_themes)
back_nut_plot <- effect_plot(back_nut_effect_sum, ci.lb, ci.ub, rowname, pal_star[1], cont_labs, back_n_labs) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  theme_black()

nat_plot <- effect_plot(nat_effect_sum, ci.lb_adjust,ci.up_adjust, duh, pal_star[2], nat_labs, nat_n_labs) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme_black()

pp_plot <- effect_plot(pp_effect_sum, ci.lb_adjust, ci.up_adjust, duh, pal_star[3], pp_labs, pp_n_labs) +
  xlab("Effect Size") +
  scale_x_continuous(limits = c(-0.5, 2),
                     labels = ~paste(., round((10^.), digits = 1), sep = "\n"),
                     name = "Effect size (log scale)\nEffect Size (response ratio)")+
  theme_black()

# Effect size plot (back nut model) ------
back_nut_plot + nat_plot + pp_plot + 
  plot_layout(ncol = 1, heights = unit(c(1.3, 6.6, 7.9), c("cm", "null"))) 

#ggsave("Output/Pub_figs/bn_effect_size_mod3.png", device = "png",
#       height = 9, width = 16, dpi = 600)


# Graph the other models! --------
# do the same thing for a supplement figure

# Lat effect size
lat_effect_sum <- coef(summary(model_full_lat)) %>%
  slice(-(3:9)) %>%
  slice(-1) %>%
  rownames_to_column() %>%
  mutate(dummy = 1, effect = estimate)

# Experimental design effect sizes
int_nat2 <- coef(summary(model_full_lat))[1,1] # extract the intercept

# edit bin labels
nat_effect_sum2 <- coef(summary(model_full_lat)) %>%
  slice(-(2)) %>%
  slice(-(6:8)) %>%
  mutate(nat_bin = factor(c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")), # add intercept to the estimated coefficient
         effect = ifelse(nat_bin != "Field: unstocked", estimate + int_nat2,
                         int_nat2), # add intercept to the estimated coefficient
         ci.lb_adjust = effect - (1.96*se), #adjust the CI based on the new estimate
         ci.up_adjust = effect + (1.96*se),
         duh = "yes") %>% 
  left_join(nat_n, by = "nat_bin")

# Primary producer effect sizes
int_pp2 <- coef(summary(model_full_lat))[1,1] # extract the intercept

# Pull rest of the coefficients, add a column for pp, and adjust the estimates to reflect the estimated effect sizes
pp_effect_sum2 <- coef(summary(model_full_lat)) %>%
  slice(-(2:6)) %>%
  mutate(pp_group = factor(c("Algae", "Coral", "Phytoplankton", "Plants")), 
         effect = ifelse(pp_group != "Algae", estimate + int_pp2,
                         int_pp2), # add intercept to the estimated coefficient
         ci.lb_adjust = effect - (1.96*se), #adjust the CI based on the new estimate
         ci.up_adjust = effect + (1.96*se),
         duh = "yes") %>% 
  left_join(pp_n, by = "pp_group")


# Create the plots
lat_labs <- c("Absolute latitude")

lat_plot <- effect_plot(lat_effect_sum, ci.lb, ci.ub, rowname, pal_star[1], lat_labs, back_n_labs) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())

nat_plot2 <- effect_plot(nat_effect_sum2, ci.lb_adjust,ci.up_adjust, duh, pal_star[2], nat_labs, nat_n_labs) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) 

pp_plot2 <- effect_plot(pp_effect_sum2, ci.lb_adjust, ci.up_adjust, duh, pal_star[3], pp_labs, pp_n_labs) +
  xlab("Effect Size") +
  scale_x_continuous(limits = c(-0.5, 2),
                     labels = ~paste(., round((10^.), digits = 1), sep = "\n"),
                     name = "Effect size (log scale)\nEffect Size (response ratio)")

# Lat effect size plot ------
lat_plot + nat_plot2 + pp_plot2 + 
  plot_layout(ncol = 1, heights = unit(c(1.3, 6.6, 7.9), c("cm", "null")))

#ggsave("Output/Pub_figs/lat_effect_size_mod2.png", device = "png",
#       height = 9, width = 16, dpi = 600)

# Temperature effect size
temp_effect_sum <- coef(summary(model_full_temp)) %>%
  slice(-(3:9)) %>%
  slice(-1) %>%
  rownames_to_column() %>%
  mutate(dummy = 1, effect = estimate)

# Experimental design effect sizes
int_nat3 <- coef(summary(model_full_temp))[1,1] # extract the intercept

# edit bin labels
nat_effect_sum3 <- coef(summary(model_full_temp)) %>%
  slice(-(2)) %>%
  slice(-(6:8)) %>%
  mutate(nat_bin = factor(c("Field: unstocked", "Field: PP stocked", "Field: animals stocked", "Field: both stocked", "Mesocosm")), # add intercept to the estimated coefficient
         effect = ifelse(nat_bin != "Field: unstocked", estimate + int_nat3,
                         int_nat3), # add intercept to the estimated coefficient
         ci.lb_adjust = effect - (1.96*se), #adjust the CI based on the new estimate
         ci.up_adjust = effect + (1.96*se),
         duh = "yes") %>% 
  left_join(nat_n, by = "nat_bin")

# Primary producer effect sizes
int_pp3 <- coef(summary(model_full_temp))[1,1] # extract the intercept

# Pull rest of the coefficients, add a column for pp, and adjust the estimates to reflect the estimated effect sizes
pp_effect_sum3 <- coef(summary(model_full_temp)) %>%
  slice(-(2:6)) %>%
  mutate(pp_group = factor(c("Algae", "Coral", "Phytoplankton", "Plants")), 
         effect = ifelse(pp_group != "Algae", estimate + int_pp3,
                         int_pp3), # add intercept to the estimated coefficient
         ci.lb_adjust = effect - (1.96*se), #adjust the CI based on the new estimate
         ci.up_adjust = effect + (1.96*se),
         duh = "yes") %>% 
  left_join(pp_n, by = "pp_group")

# Create the plots
temp_labs <- c("Temperature")

temp_plot <- effect_plot(temp_effect_sum, ci.lb, ci.ub, rowname, pal_star[1], temp_labs, back_n_labs) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())

nat_plot3 <- effect_plot(nat_effect_sum3, ci.lb_adjust,ci.up_adjust, duh, pal_star[2], nat_labs, nat_n_labs) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) 

pp_plot3 <- effect_plot(pp_effect_sum3, ci.lb_adjust, ci.up_adjust, duh, pal_star[3], pp_labs, pp_n_labs) +
  xlab("Effect Size") +
  scale_x_continuous(limits = c(-0.5, 2),
                     labels = ~paste(., round((10^.), digits = 1), sep = "\n"),
                     name = "Effect size (log scale)\nEffect Size (response ratio)")

# Temp effect size plot ------
temp_plot + nat_plot3 + pp_plot3 + 
  plot_layout(ncol = 1, heights = unit(c(1.3, 6.6, 7.9), c("cm", "null")))

#ggsave("Output/Pub_figs/temp_effect_size_mod2.png", device = "png",
#       height = 9, width = 16, dpi = 600)

# Figure 1 -----

# Figure 1a is the study site map
lat_long <- meta_data %>%
  distinct(paper_num, lat_dec, long_dec, climate) %>%
  as.data.frame() %>%
  mutate(lat = lat_dec,
         long = long_dec,
         climate = factor(climate, levels = c("Polar", "Temperate", "Subtropics", "Tropics"))) %>%
  st_as_sf(coords = c("long_dec", "lat_dec")) %>%
  st_set_crs(4326)

# Transform to a crs that uses meters as the distance
coords_jitter <- meta_coords %>%
  st_transform(3358) %>%
  st_jitter(nc_points, amount = 150000) %>%
  st_transform(4326)

# Load simple world map from package "rnaturalearth"
world <- ne_countries(scale = "medium", returnclass = "sf")

# Good sites for plotting help
# https://r-spatial.org/r/2018/10/25/ggplot2-sf.html

# https://mhallwor.github.io/_pages/Tidyverse_intro

# Define palettes
blue <- paste("#b9d1df", sep="")
pal_ew <- c("#00496f", "#0f85a0", "#ed8b00", "#dd4124")

# set seed
set.seed(987654321)

#Fig1a <- 
ggplot() +
  geom_sf(data = world, fill = "white", colour = blue) +
  geom_sf(data = coords_jitter, 
          aes(fill = climate), 
          colour = "black", pch = 21, size = 3, alpha = 0.5) +
  theme_bw() +
  theme(panel.background = element_rect(fill = blue),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 13, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 11, colour = "black")) +
  xlab("Longitude") + ylab("Latitude") + labs(fill = "Climate") +
  coord_sf(expand = FALSE)+
  scale_y_continuous(limits = c(-79, 89)) +
  scale_x_continuous(limits = c(-180, 181)) +
  scale_fill_manual(values = (pal_ew))

#ggsave("Output/Pub_figs/meta_map.png", width = 7.08661, height = 3, dpi = 600)


# Now make Fig 1b
# This is the combination of animals and pp used
set.seed(1234) # so things stop jumping around!

#Fig1b <- 
meta_data %>% 
  filter(climate != "Subtropics") %>%
  rbind(meta_data %>% filter(climate == "Subtropics")) %>%
  mutate(climate = factor(climate, 
        levels = c("Polar", "Temperate", "Subtropics", "Tropics"))) %>%
  ggplot() +
  geom_jitter(aes(c_group, pp_groupz, fill = climate, pch = nat_bin),
              size = 3.5, alpha = 0.3, width = 0.35, height = 0.35) +
  theme_classic() +
  labs(x = "Consumer taxa", y = "Primary producer taxa",
       pch = "Experiment", fill = "Climate") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 13, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 11, colour = "black")) +
  scale_fill_manual(values = (pal_ew)) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  guides(fill = guide_legend(override.aes = 
                list(size = 3.5, pch = 21, alpha = 0.4)),
         pch = guide_legend(override.aes = 
                list(size = 3.5, colour = "black", alpha = 0.9)))

#ggsave("Output/Pub_figs/taxa_region.png", device = "png",
#       height = 4, width = 7.08661, dpi = 600)

# Fig2a/Fig2b

# Plot just cont effects 
# GRAPH ALL CONT EFFECTS
cont_effect_sum <- rbind(back_nut_effect_sum, temp_effect_sum, lat_effect_sum)
# these effect tables are generated above when I'm plotting each fig by itself

# just a heinous plot
ggplot(cont_effect_sum, aes(y = dummy, x = estimate, xmin = ci.lb, xmax = ci.ub, colour = rowname)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  labs(x = "Standardized effect size", y = " ")  +
  theme_black() +
  theme(legend.position = "none") +
  scale_y_discrete(labels = cont_labs) +
  scale_colour_manual(values = rev(pal_sm)) +
  facet_wrap(~ rowname, nrow = 3, labeller = labeller(rowname = cont_labs)) + 
  theme(panel.spacing = unit(2, "lines"))
  
#ggsave("Output/Pub_figs/cont_effects.png", device = "png",
#     height = 9, width = 16, dpi = 600)
  
# Tables/summary stats ------

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

# how many papers looked at growth vs not
growth <- meta_data %>%
  group_by(growth_rate) %>%
  summarise(n = n())
  
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
  count(pp_group, pp_family, pp_genus, pp_speces) # typo in datasheet

c_used <- meta_data %>%
  count(c_group, c_family, c_genus, c_species)

c_used <- meta_data %>%
  mutate(c_activity = 
           case_when(c_group == "Bird" | c_group == "Fish" ~ "Vertebrate",
                     c_group == "Echinoderm" | c_group == "Cnidaria" ~ "Mobile Invert",
                     c_group == "Mollusc" | c_group == "Sponge"  ~ "Sessile Invert")) %>%
  count(c_activity, c_group, c_family, c_genus, c_species)

cross <- meta_data %>%
  count(pp_groupz, c_group) %>%
  mutate(percent = n/sum(n))
  
# what percentage of studies fell into each category?
# Experimental design
nat_n2 <- nat_n %>%
  mutate(percent = n/sum(n),
         delta = 0.2 - percent)

# PP group
pp_n2 <- pp_n %>%
  mutate(percent = n/sum(n))

# consumers used?
c_n <- meta_data %>%
  mutate(c_activity = case_when(c_group == "Bird" | c_group == "Fish" ~ "Vertebrate",
                                c_group == "Echinoderm" | c_group == "Cnidaria" ~ "Mobile Invert",
                                c_group == "Mollusc" | c_group == "Sponge"  ~ "Sessile Invert")) %>%
  group_by(c_activity) %>%
  summarise(n = n()) %>%
  mutate(percent = n/sum(n))   
# one mobile invert should be a sessile
# And one sessile should be mobile
# So if evens out! No prob!


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

# Which papers did I get background nut or temp proxies for? ----
bn_sum <- meta_data %>%
  group_by(nutrient_levels) %>%
  summarise(background_nut = mean(background_nut_umolL))

# Which papers did I estimate bn for?
bn_proxy <- meta_data %>%
  filter(bn_proxy == "proxy") %>%
  select(author, paper_num, nutrient_levels, year_of_study, season, bn_proxy_date, bn_proxy_dist, background_nut_umolL, comments) %>%
  distinct(paper_num, .keep_all = TRUE)
# Ok so I proxied 19 papers

temp_proxy <- meta_data %>%
  filter(temp_proxy == "proxy") %>%
  select(author, paper_num, lab_field, year_of_study, season, temp_proxy_date, temp_proxy_dist_km, temp, comments) %>%
  distinct(paper_num, .keep_all = TRUE)


# Independence and Bias -------

funnel.rma(model_all)
funnel.rma(model_full_lat)
funnel.rma(model_full_back)
funnel.rma(model_full_temp)
funnel.rma(model_null)

# A bit of a right skew

fsn(yi = yi, vi = vi, data = meta_data, type="Rosenthal") #file drawer analysis (missing studies have an effect size = 0)
# Fail-safe N: 50,340
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

sum(max_obs$n_obs) # 178
# so 10% of that would be 18 observations per paper
# so that would mean that paper # 17 and #29 have just barely too many

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
# Test for Residual Heterogeneity:
# QE(df = 171) = 2964.3669, p-val < .0001

# Test of Moderators (coefficients 2:9):
#   F(df1 = 8, df2 = 171) = 1.2423, p-val = 0.2773
# 
# Model Results:
#   
#                         estimate      se     tval   df    pval    ci.lb   ci.ub 
# intrcpt                  0.5126  0.1547   3.3137  171  0.0011   0.2072  0.8179  ** 
# lat_stand               -0.0156  0.0750  -0.2074  171  0.8359  -0.1636  0.1325     
# natural2                -0.0751  0.2264  -0.3318  171  0.7405  -0.5220  0.3718     
# natural3                 0.2667  0.2125   1.2552  171  0.2111  -0.1527  0.6862     
# natural4                 0.6122  0.3103   1.9728  171  0.0501  -0.0003  1.2248   . 
# natural5                 0.1135  0.1980   0.5733  171  0.5672  -0.2773  0.5043     
# pp_groupCoral           -0.2395  0.1971  -1.2151  171  0.2260  -0.6287  0.1496     
# pp_groupPhytoplankton   -0.3803  0.2254  -1.6867  171  0.0935  -0.8253  0.0648   . 
# pp_groupPlants          -0.3673  0.1897  -1.9363  171  0.0545  -0.7418  0.0071   .  
# 
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
# F(df1 = 8, df2 = 171) = 1.4532, p-val = 0.1778
# 
# Model Results:
#   
#                         estimate      se     tval   df    pval    ci.lb   ci.ub 
# intrcpt                  0.4542  0.1599   2.8399  171  0.0051   0.1385  0.7699  ** 
# background_nut_stand     0.0983  0.0725   1.3553  171  0.1771  -0.0449  0.2415     
# natural2                -0.0884  0.2266  -0.3899  171  0.6971  -0.5357  0.3590     
# natural3                 0.3691  0.2192   1.6833  171  0.0941  -0.0637  0.8018   . 
# natural4                 0.6605  0.3111   2.1233  171  0.0352   0.0465  1.2745   * 
# natural5                 0.1307  0.2003   0.6528  171  0.5148  -0.2646  0.5260     
# pp_groupCoral           -0.1168  0.1975  -0.5914  171  0.5551  -0.5067  0.2731     
# pp_groupPhytoplankton   -0.3669  0.2265  -1.6195  171  0.1072  -0.8140  0.0803     
# pp_groupPlants          -0.3467  0.1906  -1.8193  171  0.0706  -0.7229  0.0295   . 

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
# F(df1 = 8, df2 = 171) = 1.3100, p-val = 0.2414
# 
# Model Results:
#   
#                         estimate      se     tval   df    pval    ci.lb   ci.ub 
# intrcpt                  0.4898  0.1551   3.1577  171  0.0019   0.1836  0.7960  ** 
# temp_stand              -0.0493  0.0697  -0.7072  171  0.4804  -0.1869  0.0883     
# natural2                -0.0673  0.2259  -0.2979  171  0.7661  -0.5131  0.3785     
# natural3                 0.3069  0.2088   1.4701  171  0.1434  -0.1052  0.7190     
# natural4                 0.6428  0.3101   2.0731  171  0.0397   0.0308  1.2548   * 
# natural5                 0.1082  0.1968   0.5497  171  0.5833  -0.2803  0.4966     
# pp_groupCoral           -0.1725  0.1937  -0.8907  171  0.3743  -0.5549  0.2098     
# pp_groupPhytoplankton   -0.3836  0.2228  -1.7212  171  0.0870  -0.8234  0.0563   . 
# pp_groupPlants          -0.3699  0.1884  -1.9630  171  0.0513  -0.7419  0.0021   . 

# Dropping #17 natural2 coefficient goes from -0.2180 to 0.3455
     # Coral coefficient goes from -0.1873 to -0.4487 and p-val goes from 0.3711 to 0.0596
# Dropping #23 doesn't change anything
# Dropping #29 doesn't change anything
