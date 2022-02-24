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

source(file = "Code/theme_black_code.R")

#In the first case, object x is a fitted model object coming from the rma.uni, rma.mh, or rma.peto functions. The corresponding method is then forest.rma.


# mine data
# https://automeris.io/WebPlotDigitizer/


# Load data
# And do a million manipulations......
meta_full <- read_csv("Output/Data/nov_meta.csv") %>%
  mutate(sd_treatment = replace_na(sd_treatment, 999), # assigning a dummy to replace
         sd_control = replace_na(sd_control, 999)) %>% # assigning a dummy to replace
  mutate(sd_treatment = ifelse(sd_treatment == 999, se_treatment*sqrt(treatment_n), sd_treatment), # replacing that dummy with the estimated SD
         sd_control = ifelse(sd_control == 999, se_control*sqrt(control_n), sd_control),# replacing that dummy with the estimated SE
         natural = as.factor(natural),
         abs_lat = abs(lat_dec), # creating an absolute value lat
         climate2 = factor(ifelse(abs_lat <= 23.45, "Tropics",
                           ifelse(abs_lat <= 40, "Subtropics",
                                  ifelse(abs_lat <= 60, "Temperate",
                                         "Polar"))),
         levels = c("Tropics", "Subtropics", "Temperate", "Polar")), # New column where I don't distingiush between North subtropics/temperates and South subtropics/temperates
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
# Another paper set 23.5 - 23.5, and set temperete as 30-40

# Calculate effect sizes -------

##### I think I need to manually exclude the negative means because it's still calculating effect sizes for those
meta_no_neg <- meta_full %>%
  filter(mean_control > 0) %>%
  filter(mean_treatment > 0)

# We lose 5 rows when we remove 0's and 7 when we remove negatives (12 total)
# should I switch to Hedges g????
# Should I add like 0.0000000000001 to all the 0's?
 # NO DO NOT DO THIS IT RUINS EVERYTHING EVERY STATISTICAL ANALYSIS BECOMES P = 0.99

# Try calculating Hedges g 
meta_g <- meta_full %>%
  escalc(measure = "SMD", data = meta_no_neg,
         m1i = mean_treatment, sd1i = sd_treatment, 
         m2i = mean_control, sd2i = sd_control,
         n1i = treatment_n, n2i = control_n) %>%
  transmute(indiv = indiv,
            effect_hedges_g = yi,
            var_hedges_g = vi)

# OK so when I use Hedges g for some reason chl a goes not significant. 
# I think we should stick with ROM??? That's the effect size we initially chose

meta <- meta_no_neg %>%
  escalc(measure = "ROM", data = meta_no_neg,
       m1i = mean_treatment, sd1i = sd_treatment, 
       m2i = mean_control, sd2i = sd_control,
       n1i = treatment_n, n2i = control_n) %>%
# ROM for log transformed ratio of means
# log(treatment/control)
  left_join(meta_g, by = "indiv")

# Compare hedges g vs ROM
ggplot(meta, aes(yi, effect_hedges_g)) +
  geom_point()

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
  select(author, title, year, lat_dec, long_dec, climate_zone, climate2)

meta_all_papers <- meta_reduced %>%
  select(author, title, year, lat_dec, long_dec, climate_zone, climate2) %>%
  rbind(meta_rejected) # df with all the papers I'm using and the papers I can't use but that still mostly meet the criteria

# Add chla data to orig meta dataframe  -------
# this file was writted in the Chla_extraction.R code file
chla_all_values <- read_csv("Output/Data/all_extracted_chla.csv")


# What's the max and min chla value?
max <- max(chla_all_values$chla_est) # 26.73629
min<- min(chla_all_values$chla_est) # 0.05139377

# How big should each bin be?
bin_size <- max/5
upper_bin1 <- bin_size
upper_bin2 <- bin_size*2
upper_bin3 <- bin_size*3
upper_bin4 <- bin_size*4
upper_bin5 <- bin_size*5

# use this so that it still works even if the max chla value changes

# add chl a estimates to df
# add a binned chl a column
meta_final <- meta %>%
  left_join(chla_all_values, by = "indiv") %>%
  mutate(chla_bin = ifelse(chla_est <= upper_bin1, "1",
                      ifelse(chla_est <= upper_bin2, "2",
                             ifelse(chla_est <= upper_bin3, "3", 
                                 ifelse(chla_est <= upper_bin4, "4", "5")))))

ggplot(meta_final, aes(chla_bin, chla_est, colour = chla_bin)) +
  geom_jitter() +
  theme_classic()

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


#ggsave("Output/Figures/climate_zone_bar", device = "png",
#      height = 9, width = 16, dpi = 400)


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

# chla vs effect size?
ggplot(meta_final, aes(chla_est, yi)) +
  geom_point() +
  geom_smooth(method = lm)

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
  

# stats???------

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


# Multilevel models ------
# There is a special way through which we can tell rma.mv to assume such nested random effects. We do this using a slash (/) to separate the higher- and lower-level grouping variable. To the left of /, we put in the level 3 (cluster) variable. To the right, we insert the lower-order variable nested in the larger cluster. Therefore, the general structure of the formula looks like this: ~ 1 | cluster/effects_within_cluster

# Lets look at latitude first ----
model_lat_paper <- rma.mv(yi = yi,
              V = vi,
              slab = author,
              data = meta_final,
              method = "REML",
              mods = ~ climate_zone,
              random = ~ 1 | paper_num/indiv, 
              test = "t")

summary(model_lat_paper)
# Look at Test of Moderators to detemine whether there are significant differences between the subgroups
# p-value is 0.9101 so NO, there is NO difference!
# The Model Results are printed within a meta-regression framework. This means that we cannot directly extract the estimates in order to obtain the pooled effect sizes within subgroups.

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

# Try using the 4 climate zones instead ----
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
n_labs_cli <- c("n = 28", "n = 70", "n = 74", "n = 2")
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


# Now try with chl a -----
model_chla_paper <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_final,
                          method = "REML",
                          mods = ~ chla_bin,
                          random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                          test = "t")

summary(model_chla_paper)
# Look at Test of Moderators to detemine whether there are significant differences between the subgroups
# p-value is 0.0148 so yes, there is a difference!
# The Model Results are printed within a meta-regression framework. This means that we cannot directly extract the estimates in order to obtain the pooled effect sizes within subgroups.

# How many papers in each bin?
chla_n <- meta_final %>%
  group_by(chla_bin) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:5,
         chla_bin = factor(c("Very low", "Low", "Medium", "High", "Very high")))

# Pull coefficient estimates so I can graph them
intercept_chla <- coef(summary(model_chla_paper))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes
chla_effect_sum <- coef(summary(model_chla_paper)) %>%
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
n_labs_chla <- c("n = 111", "n = 30", "n = 7", "n = 6", "n = 20")

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

# OK OK OK SO LOWER CHL A HAS A HIGHER EFFECT!!!!! WOOOOOOOW


# do a permutation test to make sure we haven't overfit the model 
#permutest(model_chla_paper)

# How about natural-ness of experiment -----
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

nat_n_labs <- c("n = 72", "n = 14", "n = 16", "n = 10", "n = 62")

ggplot(nat_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = nat_bin)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size (ROM)", y = " ") +
  scale_y_continuous(breaks = 1:length(nat_labs),
                     labels = nat_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(nat_n_labs),
                                         labels = nat_n_labs)) +
  scale_colour_manual(values = rev(pal_star)) +
  theme_black()+ 
  theme(legend.position = "none")

#ggsave("Output/Figures/effects_nat", device = "png",
#     height = 9, width = 16, dpi = 400)

# How about consumer activity level -----
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

act_n_labs <- c("n = 46", "n = 21", "n = 107")

ggplot(act_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = c_activity)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size (ROM)", y = " ") +
  scale_y_continuous(breaks = 1:length(act_labs),
                     labels = act_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(act_n_labs),
                                         labels = act_n_labs)) +
  scale_colour_manual(values = rev(pal_star)) +
  theme_black()+ 
  theme(legend.position = "none")

#ggsave("Output/Figures/effects_act", device = "png",
#     height = 9, width = 16, dpi = 400)


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


# Independence and Bias -------

funnel.rma(model_lat_paper)
funnel.rma(model_chla_paper)
funnel.rma(model_nat_paper)
# A bit of a right skew

fsn(yi = yi, vi = vi, data = meta_final, type="Rosenthal") #file drawer analysis (missing studies have an effect size = 0)
# Fail-safe N: 35,685
# That seems like quite a large number?????

# Check for influential studies? ------

# make a df where each row is a paper, and num_obs is the number of observations that paper contributed
max_obs <- meta_final %>%
  group_by(paper_num) %>%
  dplyr::summarise(num_obs = as.numeric(n()))

hist(max_obs$num_obs)

# there's 174 data points
# so 10% of that would be 17.4 observations per paper
# so that would mean that paper # 17 and #23 and #31 have too many. try running things without those two

# Ok so let's remove those papers and rerun the stats???

meta_f <- meta_final %>%
  filter(paper_num != 23)

# What's the max and min chla value?
max <- max(meta_f$chla_est) # 21.48487
min<- min(meta_f$chla_est) # 0.05139377

# How big should each bin be?
bin_size <- max/5
upper_bin1 <- bin_size
upper_bin2 <- bin_size*2
upper_bin3 <- bin_size*3
upper_bin4 <- bin_size*4
upper_bin5 <- bin_size*5

meta_f <- meta_f %>%
  mutate(chla_bin = ifelse(chla_est <= upper_bin1, "1",
                           ifelse(chla_est <= upper_bin2, "2",
                                  ifelse(chla_est <= upper_bin3, "3", 
                                         ifelse(chla_est <= upper_bin4, "4", "5")))))


# Rerun lat ----
re_model_lat_paper <- rma.mv(yi = yi,
                          V = vi,
                          slab = author,
                          data = meta_f,
                          method = "REML",
                          mods = ~ climate_zone,
                          random = ~ 1 | paper_num/indiv, 
                          test = "t")

summary(re_model_lat_paper)
# Removed # 17, p-value is 0.9187 versus 0.9101
# Removed # 23, p-value is 0.8921 versus 0.9101
# Removed # 23, p-value is 0.7748 versus 0.9101


# # Pull coefficient estimates so I can graph them
re_intercept_lat <- coef(summary(re_model_lat_paper))[1,1] # extract the intercept

# How many papers in each bin?
re_lat_n <- meta_f %>%
  group_by(climate_zone) %>%
  summarise(n = n()) %>%
  rename(climate = climate_zone) %>%
  mutate(dummy = 1:6)

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes
re_lat_effect_sum <- coef(summary(re_model_lat_paper)) %>%
  mutate(climate = factor(c("ST", "SST", "T", "NST", "NT", "NP"),
                          levels = c("ST", "SST", "T", "NST", "NT", "NP")),
         effect = ifelse(climate != "ST", estimate + re_intercept_lat,
                         re_intercept_lat),
         ci.lb = ifelse(climate != "ST", ci.lb + re_intercept_lat,
                        ci.lb),
         ci.ub = ifelse(climate != "ST", ci.ub + re_intercept_lat,
                        ci.ub)) %>%# add intercept to other coefficients
  left_join(re_lat_n, by = "climate")

# GRAPH THIS
re_n_labs_lat <- c("n = 5", "n = 10", "n = 28", "n = 33", "n = 42", "n = 2")

ggplot(re_lat_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = climate)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept=0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size (ROM)", y = "Climate zone") +
  scale_y_continuous(breaks = 1:length(climate_labs),
                     labels = climate_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(re_n_labs_lat),
                                         labels = re_n_labs_lat)) +
  scale_colour_manual(values = rev(pal))  +
  theme_black()+ 
  theme(legend.position = "none")

# Removed # 17 all the bars still overlap O
# Removed # 23 all the bars still overlap O
# Removed # 31 all the bars still overlap O


# Rerun with 4 climate zones  ----
re_model_climate_paper <- rma.mv(yi = yi,
                              V = vi,
                              slab = author,
                              data = meta_f,
                              method = "REML",
                              mods = ~ climate2,
                              random = ~ 1 | paper_num/indiv, 
                              test = "t")

summary(re_model_climate_paper)
# Removed # 17, p-value is 0.7605 versus 0.6538
# Removed # 23, p-value is 0.7563 versus 0.6538
# Removed # 31, p-value is 0.1750 versus 0.6538


# # Pull coefficient estimates so I can graph them
re_intercept_cli <- coef(summary(re_model_climate_paper))[1,1] # extract the intercept

# How many papers in each bin?
re_cli_n <- meta_f %>%
  group_by(climate2) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:4)

# Now pull the rest of the coefficients, add a column for climate, and adjust the estimates to reflect the estimated effect sizes
re_cli_effect_sum <- coef(summary(re_model_climate_paper)) %>%
  mutate(climate2 = factor(c("Tropics", "Subtropics", "Temperate", "Polar"),
                           levels = c("Tropics", "Subtropics", "Temperate", "Polar")),
         effect = ifelse(climate2 != "Tropics", estimate + re_intercept_cli,
                         re_intercept_cli),
         ci.lb = ifelse(climate2 != "Tropics", ci.lb + re_intercept_cli,
                        ci.lb),
         ci.ub = ifelse(climate2 != "Tropics", ci.ub + re_intercept_cli,
                        ci.ub)) %>% # add intercept to other coefficients
  left_join(re_cli_n, by = "climate2")


# GRAPH THIS
n_labs_cli <- c("n = 28", "n = 70", "n = 56", "n = 2")
cli_labs <- c("Tropics", "Subtropics", "Temperate", "Polar")

ggplot(re_cli_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = climate2)) +
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

# Tropics and subtropics still positive, temp still intercepts 0


# Rerun chl a -----
re_model_chla_paper <- rma.mv(yi = yi,
                           V = vi,
                           slab = author,
                           data = meta_f,
                           method = "REML",
                           mods = ~ chla_bin,
                           random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
                           test = "t")

summary(re_model_chla_paper)
# Remove # 17 p-value is 0.0090 versus 0.0148 
# Remove # 23 p-value is 0.0633 versus 0.0148 
# Remove # 23 and don't change bins 0.0818 versus 0.0148 
# Remove # 31 p-value is 0.0491 versus 0.0148 



# How many papers in each bin?
re_chla_n <- meta_f %>%
  group_by(chla_bin) %>%
  summarise(n = n()) %>%
  mutate(dummy = 1:5,
         chla_bin = factor(c("Very low", "Low", "Medium", "High", "Very high")))

# Pull coefficient estimates so I can graph them
re_intercept_chla <- coef(summary(re_model_chla_paper))[1,1] # extract the intercept

# Now pull the rest of the coefficients, add a column for chla, and adjust the estimates to reflect the estimated effect sizes
re_chla_effect_sum <- coef(summary(re_model_chla_paper)) %>%
  mutate(chla_bin = factor(c("Very low", "Low", "Medium", "High", "Very high"),
                           levels = c("Very low", "Low", "Medium", "High", "Very high")),
         effect = ifelse(chla_bin != "Very low", estimate + re_intercept_chla,
                         re_intercept_chla),
         ci.lb = ifelse(chla_bin != "Very low", ci.lb + re_intercept_chla,
                        ci.lb),
         ci.ub = ifelse(chla_bin != "Very low", ci.ub + re_intercept_chla,
                        ci.ub)) %>%# add intercept to other coefficients
  left_join(re_chla_n, by = "chla_bin")

# GRAPH THIS?
chla_bin_labs <- c("Very low", "Low", "Medium", "High", "Very high")
re_n_labs_chla <- c("n = 111", "n = 12", "n = 7", "n = 6", "n = 20")

ggplot(re_chla_effect_sum, aes(y = dummy, x = effect, xmin = ci.lb, xmax = ci.ub, colour = chla_bin)) +
  geom_point(size = 6) +
  geom_errorbar(width = .5, size = 2) +
  geom_vline(xintercept = 0, color="white", linetype="dashed") +
  theme_classic() +
  labs(x = "Effect size (ROM)", y = "Chlorophyll a") +
  scale_y_continuous(breaks = 1:length(chla_bin_labs),
                     labels = chla_bin_labs,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(re_n_labs_chla),
                                         labels = re_n_labs_chla)) +
  scale_color_brewer(palette = 'Greens') +
  theme_black()+ 
  theme(legend.position = "none")

# Removed # 17 looks basically the same
# Removed # 23 now the medium is positive, and the bins are different (bc removing 23 changes the bins bc it was the max chl a )
# Removed # 23 didn't change bin sizes and the graph basically looks the same
# Removed # 31 looks basically the same


# When I remove all three chl a increases with chl a until we get to the high chl a and then it's a "bad" thing...
# But pp in areas with lots of nutrients are theorized to have better uptake abilities
# so maybe if there are medium-ish nutrients around, then the pp will have better uptake pathways so they can better take advantage of the extra pee in the water
# 


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

