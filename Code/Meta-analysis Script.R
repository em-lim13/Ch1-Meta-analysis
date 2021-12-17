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

# packages for raster stuff?
#library(ncdf4)
#library(RNetCDF)
#library(tidync)

#library(ncmeta)
#library(stars)
#library(tidyr)
library(raster)
library(rgdal)
library(sf)
library(sp)


#In the first case, object x is a fitted model object coming from the rma.uni, rma.mh, or rma.peto functions. The corresponding method is then forest.rma.


# mine data
# https://automeris.io/WebPlotDigitizer/


# Load data
meta_full <- read_csv("Output/Data/nov_meta.csv") %>%
  mutate(climate_zone = factor(climate_zone, levels = c("SP", "ST", "SST", "T", "NST", "NT", "NP"))) %>%
  mutate(sd_treatment = replace_na(sd_treatment, 999),
         sd_control = replace_na(sd_control, 999)) %>%
  mutate(sd_treatment = ifelse(sd_treatment == 999, se_treatment*sqrt(treatment_n), sd_treatment),
         sd_control = ifelse(sd_control == 999, se_control*sqrt(control_n), sd_control))


meta_no_neg <- meta_full %>%
  filter(mean_control >= 0) %>%
  filter(mean_treatment >= 0)

meta <- meta_no_neg %>%
  escalc(measure = "ROM", data = meta_no_neg,
       m1i = mean_treatment, sd1i = sd_treatment, 
       m2i = mean_control, sd2i = sd_control,
       n1i = treatment_n, n2i = control_n) 

##### I think I need to manually exclude the negative means because it's still calculating effect sizes for those

# Make a data frame where there's only one row for each paper to graph where the papers are from
meta_reduced <- meta_full %>%
  filter(obs_num == 1) 


# Load papers that I can't use but that mostly pass the criteria
meta_rejected <- read_csv("Output/Data/rejected_papers.csv")

meta_all_papers <- rbind(meta_reduced, meta_rejected) # df with all the papers I'm using and the papers I can't use but that still mostly meet the criteria

# Graphing --------
# Make palette
pal <- pnw_palette("Bay", 6)

#using meta_reduced so each paper is only counted once
ggplot(meta_all_papers, aes(y = climate_zone, fill = climate_zone)) +
  geom_bar(stat="count") +
  scale_fill_manual(values = rev(pal)) +
  scale_y_discrete(labels = c("South Temperate","South Subtropics", "Tropics", "North Subtropics", "North Temperate", "North Polar")) +
  labs( y = "Climate Zone", x = "Number of Papers") + 
  theme_black() + theme(legend.position = "none")


#ggsave("Output/Figures/climate_zone_bar", device = "png",
#       height = 9, width = 16, dpi = 400)


# Look at effect sizes????
ggplot(meta, aes(x = yi, y = lat_dec,  colour = climate_zone)) +
  geom_point(size = 3) +
  scale_colour_manual(values = rev(pal), labels = c("South Temperate","South Subtropics", "Tropics", "North Subtropics", "North Temperate", "North Polar"), guide = guide_legend(reverse = TRUE)) +
  labs( y = "Lattitude", x = "Effect size", colour = "Climate Zone") + 
  geom_vline(xintercept = 0, linetype="dotted", colour = "white", size = 1.6) +
  theme_black()


#ggsave("Output/Figures/response_lat", device = "png",
#       height = 9, width = 16, dpi = 400)


ggplot(meta, aes(x = yi, y = lat_dec,  colour = climate_zone)) +
  geom_point(size = 3) +
  scale_colour_manual(values = rev(pal), guide = FALSE) +
  labs( y = "Lattitude", x = "Effect size", colour = "Climate Zone") + 
  geom_vline(xintercept = 0, linetype="dotted", colour = "white", size = 1.6) +
  theme_black()

#ggsave("Output/Figures/response_lat_noleg", device = "png",
#      height = 9, width = 16, dpi = 400)

# stats???------

# Full versus reduced model comparisons in terms of fit statistics and likelihoods can be obtained with anova.rma
# Permutation tests for the model coefficient(s) can be obtained with permutest.rma.uni. Tests and confidence intervals based on (cluster) robust methods can be obtained with robust.rma.uni.
# Tests for funnel plot asymmetry (which may be indicative of publication bias) can be obtained with ranktest.rma and regtest.rma. For models without moderators, the trimfill.rma.uni method can be used to carry out a trim and fill analysis and hc.rma.uni provides a random-effects model analysis that is more robust to publication bias (based on the method by Henmi & Copas, 2010).

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


# Try using tutorial? -----
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



# Try using "Doing Meta-Analysis in R" -------
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/metareg.html

# If you don't want to compare the models directly using the anova function, use the "REML" (restricted maximum likelihood) 
model_control <- rma(yi = yi,
                 vi = vi,
                 data = meta,
                 method = "ML",
                 test = "knha")

model_lat <- rma(yi = yi,
                       vi = vi,
                       data = meta,
                       method = "ML",
                       mods = ~ lat_dec,
                       test = "knha")

# is the model with lat better than the null model? NO
anova(model_control, model_lat)

# Multilevel model
# There is a special way through which we can tell rma.mv to assume such nested random effects. We do this using a slash (/) to separate the higher- and lower-level grouping variable. To the left of /, we put in the level 3 (cluster) variable. To the right, we insert the lower-order variable nested in the larger cluster. Therefore, the general structure of the formula looks like this: ~ 1 | cluster/effects_within_cluster

model_lat_paper <- rma.mv(yi = yi,
              V = vi,
              slab = author,
              data = meta,
              method = "REML",
              mods = ~ lat_dec,
              random = ~ 1 | paper_num/indiv, # Need to decide if I need to nest or not
              test = "t")


summary(model_lat_paper)

# need to convert z estimate to a 'normal correlation
convert_z2r(-0.0010)

# do a permutation test to make sure we haven't overfit the model 
#permutest(model_lat_paper)


# Chl a data ------

# Seasonal chl a data
# Download chl a data from here
# https://oceancolor.gsfc.nasa.gov/l3/order/

# group meta data by season and year so I know which years to download
meta_years <- meta %>%
  group_by(year_of_study, season) %>%
  summarise(n = n())

# basically want a csv for each time point with lat long for each study in that time point

# Define coordinate reference system
prj4string <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
my.projection <- st_crs(prj4string)

# SUMMER 2004 -----

# Load a chl a file
# chl a value is in mg/m^2
raster_summer2004 <- raster("Data/Rasters/A20041732004264.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2004 <- meta %>%
  filter(season == "spring" & year_of_study == "2005") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer2004 <- meta_summer2004 %>%
  select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2004 <- raster::extract(raster_summer2004, sf_summer2004, df = TRUE) %>%
  bind_cols(meta_summer2004) %>%
  select(-ID)

# SPRING 2005 -----

# Load chl a
raster_spring2005 <- raster("Data/Rasters/A20050802005171.L3m_SNSP_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_spring2005 <- meta %>%
  filter(season == "spring" & year_of_study == "2005") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_spring2005 <- meta_spring2005 %>%
  select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_spring2005 <- raster::extract(raster_spring2005, sf_spring2005, df = TRUE) %>%
  bind_cols(meta_spring2005) %>%
  select(-ID)

# SUMMER 2005 -------

# Load chl a
raster_summer2005 <- raster("Data/Rasters/A20051722005263.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2005 <- meta %>%
  filter(season == "summer" & year_of_study == "2005") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer2005 <- meta_summer2005 %>%
  select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2005 <- raster::extract(raster_summer2005, sf_summer2005, df = TRUE) %>%
  bind_cols(meta_summer2005) %>%
  select(-ID)


# Old way to get shapefiles ------

# Load chl a
setwd("~/Downloads")
raster_summer2004 <- raster("Data/Rasters/A20041732004264.L3m_SNSU_CHL_chlor_a_4km.nc")
setwd("~/Documents/MSc/Ch2-Meta-analysis")

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
