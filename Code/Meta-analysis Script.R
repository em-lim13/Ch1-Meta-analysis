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
  theme_black() +
  ylim(c(-75, 75))


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

# Couldn't find dates for Uthicke 2001, Cognie and Barille 1999, and Arzul et al 2001

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
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2004 <- raster::extract(raster_summer2004, sf_summer2004, df = TRUE) %>%
  bind_cols(meta_summer2004) %>%
  dplyr::select(-ID)

# SPRING 2005 -----
# Load chl a
raster_spring2005 <- raster("Data/Rasters/A20050802005171.L3m_SNSP_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_spring2005 <- meta %>%
  filter(season == "spring" & year_of_study == "2005") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_spring2005 <- meta_spring2005 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_spring2005 <- raster::extract(raster_spring2005, sf_spring2005, df = TRUE) %>%
  bind_cols(meta_spring2005) %>%
  dplyr::select(-ID)

# SUMMER 2005 -------
# Load chl a
raster_summer2005 <- raster("Data/Rasters/A20051722005263.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2005 <- meta %>%
  filter(season == "summer" & year_of_study == "2005") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer2005 <- meta_summer2005 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2005 <- raster::extract(raster_summer2005, sf_summer2005, df = TRUE) %>%
  bind_cols(meta_summer2005) %>%
  dplyr::select(-ID)

# WINTER 2005 -------
# Load chl a
raster_winter2005 <- raster("Data/Rasters/A20053552006079.L3m_SNWI_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_winter2005 <- meta %>%
  filter(season == "winter" & year_of_study == "2005") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_winter2005 <- meta_winter2005 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_winter2005 <- raster::extract(raster_winter2005, sf_winter2005, df = TRUE) %>%
  bind_cols(meta_winter2005) %>%
  dplyr::select(-ID)

# SUMMER 2006 -------
# Load chl a
raster_summer2006 <- raster("Data/Rasters/A20061722006263.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2006 <- meta %>%
  filter(season == "summer" & year_of_study == "2006") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer2006 <- meta_summer2006 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2006 <- raster::extract(raster_summer2006, sf_summer2006, df = TRUE) %>%
  bind_cols(meta_summer2006) %>%
  dplyr::select(-ID)

# SUMMER 2009 -------
# Load chl a
raster_summer2009 <- raster("Data/Rasters/A20091722009263.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2009 <- meta %>%
  filter(season == "summer" & year_of_study == "2009") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer2009 <- meta_summer2009 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2009 <- raster::extract(raster_summer2009, sf_summer2009, df = TRUE) %>%
  bind_cols(meta_summer2009) %>%
  dplyr::select(-ID)

# SUMMER 2010 -------
# Load chl a
raster_summer2010 <- raster("Data/Rasters/A20101722010263.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2010 <- meta %>%
  filter(season == "summer" & year_of_study == "2010") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer2010 <- meta_summer2010 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2010 <- raster::extract(raster_summer2010, sf_summer2010, df = TRUE) %>%
  bind_cols(meta_summer2010) %>%
  dplyr::select(-ID)

# WINTER 2010 -------
# Load chl a
raster_winter2010 <- raster("Data/Rasters/A20103552011079.L3m_SNWI_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_winter2010 <- meta %>%
  filter(season == "winter" & year_of_study == "2010") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_winter2010 <- meta_winter2010 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_winter2010 <- raster::extract(raster_winter2010, sf_winter2010, df = TRUE) %>%
  bind_cols(meta_winter2010) %>%
  dplyr::select(-ID)

# SUMMER 2011 -------
# Load chl a
raster_summer2011 <- raster("Data/Rasters/A20111722011263.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2011 <- meta %>%
  filter(season == "summer" & year_of_study == "2011") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer2011 <- meta_summer2011 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2011 <- raster::extract(raster_summer2011, sf_summer2011, df = TRUE) %>%
  bind_cols(meta_summer2011) %>%
  dplyr::select(-ID)

# FALL 2012 -------
# Load chl a
raster_fall2012 <- raster("Data/Rasters/A20122652012355.L3m_SNAU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_fall2012 <- meta %>%
  filter(season == "fall" & year_of_study == "2012") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_fall2012 <- meta_fall2012 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_fall2012 <- raster::extract(raster_fall2012, sf_fall2012, df = TRUE) %>%
  bind_cols(meta_fall2012) %>%
  dplyr::select(-ID)

# SUMMER 2012 -------
# Load chl a
raster_summer2012 <- raster("Data/Rasters/A20121732012264.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2012 <- meta %>%
  filter(season == "summer" & year_of_study == "2012") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer2012 <- meta_summer2012 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2012 <- raster::extract(raster_summer2012, sf_summer2012, df = TRUE) %>%
  bind_cols(meta_summer2012) %>%
  dplyr::select(-ID)

# SUMMER 2013 -------
# Load chl a
raster_summer2013 <- raster("Data/Rasters/A20131722013263.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2013 <- meta %>%
  filter(season == "summer" & year_of_study == "2013") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer2013 <- meta_summer2013 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2013 <- raster::extract(raster_summer2013, sf_summer2013, df = TRUE) %>%
  bind_cols(meta_summer2013) %>%
  dplyr::select(-ID)

# WINTER 2013 -------
# Load chl a
raster_winter2013 <- raster("Data/Rasters/A20133552014079.L3m_SNWI_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_winter2013 <- meta %>%
  filter(season == "winter" & year_of_study == "2013") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_winter2013 <- meta_winter2013 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_winter2013 <- raster::extract(raster_winter2013, sf_winter2013, df = TRUE) %>%
  bind_cols(meta_winter2013) %>%
  dplyr::select(-ID)

# SPRING 2015 -------
# Load chl a
raster_spring2015 <- raster("Data/Rasters/A20150802015171.L3m_SNSP_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_spring2015 <- meta %>%
  filter(season == "spring" & year_of_study == "2015") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_spring2015 <- meta_spring2015 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_spring2015 <- raster::extract(raster_spring2015, sf_spring2015, df = TRUE) %>%
  bind_cols(meta_spring2015) %>%
  dplyr::select(-ID)

# SUMMER 2015 -------
# Load chl a
raster_summer2015 <- raster("Data/Rasters/A20151722015263.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2015 <- meta %>%
  filter(season == "summer" & year_of_study == "2015") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer2015 <- meta_summer2015 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2015 <- raster::extract(raster_summer2015, sf_summer2015, df = TRUE) %>%
  bind_cols(meta_summer2015) %>%
  dplyr::select(-ID)

# SPRING 2016 -------
# Load chl a
raster_spring2016 <- raster("Data/Rasters/A20160812016172.L3m_SNSP_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_spring2016 <- meta %>%
  filter(season == "spring" & year_of_study == "2016") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_spring2016 <- meta_spring2016 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_spring2016 <- raster::extract(raster_spring2016, sf_spring2016, df = TRUE) %>%
  bind_cols(meta_spring2016) %>%
  dplyr::select(-ID)

# WINTER 2017 -------
# Load chl a
raster_winter2017 <- raster("Data/Rasters/A20173552018079.L3m_SNWI_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_winter2017 <- meta %>%
  filter(season == "winter" & year_of_study == "2017") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_winter2017 <- meta_winter2017 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_winter2017 <- raster::extract(raster_winter2017, sf_winter2017, df = TRUE) %>%
  bind_cols(meta_winter2017) %>%
  dplyr::select(-ID)

# FALL 2019 -------

# Load chl a
raster_fall2019 <- raster("Data/Rasters/A20192642019354.L3m_SNAU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_fall2019 <- meta %>%
  filter(season == "fall" & year_of_study == "2019") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_fall2019 <- meta_fall2019 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_fall2019 <- raster::extract(raster_fall2019, sf_fall2019, df = TRUE) %>%
  bind_cols(meta_fall2019) %>%
  dplyr::select(-ID)

# 2000 - 2005 -------
# Load chl a from 2002, 2003, 2004, and 2005 (project started in 2002)
raster_2002 <- raster("Data/Rasters/A20020012002365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

raster_2003 <- raster("Data/Rasters/A20030012003365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

raster_2004 <- raster("Data/Rasters/A20040012004366.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

raster_2005 <- raster("Data/Rasters/A20050012005365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_2000_2005 <- meta %>%
  filter(year_of_study == "2000-2005") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_2000_2005 <- meta_2000_2005 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the chl a data
temp_2002a <- raster::extract(raster_2002, sf_2000_2005, df = TRUE) %>%
  rename(chla_2002 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_2003 <- raster::extract(raster_2003, sf_2000_2005, df = TRUE) %>%
  rename(chla_2003 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_2004 <- raster::extract(raster_2004, sf_2000_2005, df = TRUE) %>%
  rename(chla_2004 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_2005 <- raster::extract(raster_2005, sf_2000_2005, df = TRUE) %>%
  rename(chla_2005 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

# now bind all these yearly averages into the same dataframe and average across all years
chla_2000_2005 <- bind_cols(temp_2002a, temp_2003, temp_2004, temp_2005) %>%
  mutate(Chlorophyll.Concentration..OCI.Algorithm = rowMeans(cbind(chla_2002, chla_2003, chla_2004, chla_2005))) %>%
  dplyr::select(Chlorophyll.Concentration..OCI.Algorithm) %>%
  bind_cols(meta_2000_2005)

# 2001 - 2002 -------
# Load chl a (just from 2002 bc that's when the project started)
raster_2002 <- raster("Data/Rasters/A20020012002365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_2002 <- meta %>%
  filter(year_of_study == "2001-2002") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_2002 <- meta_2002 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_2002 <- raster::extract(raster_2002, sf_2002, df = TRUE) %>%
  bind_cols(meta_2002) %>%
  dplyr::select(-ID)

# 2003 - 2005 -------
# Load chl a from 2003, 2004, and 2005
raster_2003 <- raster("Data/Rasters/A20030012003365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

raster_2004 <- raster("Data/Rasters/A20040012004366.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

raster_2005 <- raster("Data/Rasters/A20050012005365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_2003_2005 <- meta %>%
  filter(year_of_study == "2003-2005") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_2003_2005 <- meta_2003_2005 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the chl a data
temp_2003b <- raster::extract(raster_2003, sf_2003_2005, df = TRUE) %>%
  rename(chla_2003 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_2004b <- raster::extract(raster_2004, sf_2003_2005, df = TRUE) %>%
  rename(chla_2004 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_2005b <- raster::extract(raster_2005, sf_2003_2005, df = TRUE) %>%
  rename(chla_2005 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

# now bind all these yearly averages into the same dataframe and average across all years
chla_2003_2005 <- bind_cols(temp_2003b, temp_2004b, temp_2005b) %>%
  mutate(Chlorophyll.Concentration..OCI.Algorithm = rowMeans(cbind(chla_2003, chla_2004, chla_2005))) %>%
  dplyr::select(Chlorophyll.Concentration..OCI.Algorithm) %>%
  bind_cols(meta_2003_2005)

# 2008 --------
# Load chl a (yearly average for 2008)
raster_2008 <- raster("Data/Rasters/A20080012008366.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one year
# CAREFUL BC IF THERE'S ANOTHER STUDY IN 2008 THAT DOES HAVE A SPECIFIC SEASON THIS WOULD CATCH THAT ONE TOO
meta_2008 <- meta %>%
  filter(year_of_study == "2008") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_2008 <- meta_2008 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_2008 <- raster::extract(raster_2008, sf_2008, df = TRUE) %>%
  bind_cols(meta_2008) %>%
  dplyr::select(-ID)

# 2013 - 2014 -------
# Load yearly average chl a from 2013, 2014
raster_2013 <- raster("Data/Rasters/A20130012013365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

raster_2014 <- raster("Data/Rasters/A20140012014365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_2013_2014 <- meta %>%
  filter(year_of_study == "2013-2014") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_2013_2014 <- meta_2013_2014 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the chl a data
temp_2013 <- raster::extract(raster_2013, sf_2013_2014, df = TRUE) %>%
  rename(chla_2013 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_2014 <- raster::extract(raster_2014, sf_2013_2014, df = TRUE) %>%
  rename(chla_2014 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

# now bind all these yearly averages into the same dataframe and average across all years
chla_2013_2014 <- bind_cols(temp_2013, temp_2014) %>%
  mutate(Chlorophyll.Concentration..OCI.Algorithm = rowMeans(cbind(chla_2013, chla_2014))) %>%
  dplyr::select(Chlorophyll.Concentration..OCI.Algorithm) %>%
  bind_cols(meta_2013_2014)

# 2015 - 2016 -------
# Load yearly average chl a from 2015, 2016
raster_2015 <- raster("Data/Rasters/A20150012015365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

raster_2016 <- raster("Data/Rasters/A20160012016366.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_2015_2016 <- meta %>%
  filter(year_of_study == "2015-2016") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_2015_2016 <- meta_2015_2016 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the chl a data
temp_2015 <- raster::extract(raster_2015, sf_2015_2016, df = TRUE) %>%
  rename(chla_2015 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_2016 <- raster::extract(raster_2016, sf_2015_2016, df = TRUE) %>%
  rename(chla_2016 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

# now bind all these yearly averages into the same dataframe and average across all years
chla_2015_2016 <- bind_cols(temp_2015, temp_2016) %>%
  mutate(Chlorophyll.Concentration..OCI.Algorithm = rowMeans(cbind(chla_2015, chla_2016))) %>%
  dplyr::select(Chlorophyll.Concentration..OCI.Algorithm) %>%
  bind_cols(meta_2015_2016)

# LEFTOVER PAPERS -------

# Stick all the completed chla dataframes together????
chla_prelim <- rbind(chla_2000_2005, chla_2002, chla_2003_2005, chla_2008, chla_2013_2014, chla_2015_2016, chla_fall2012, chla_fall2019, chla_spring2005, chla_spring2015, chla_spring2016, chla_summer2004, chla_summer2005, chla_summer2006, chla_summer2009, chla_summer2010, chla_summer2011, chla_summer2012, chla_summer2013, chla_summer2015, chla_winter2005, chla_winter2010, chla_winter2013, chla_winter2017) %>% #24 dfs 
  rename(chla_est = Chlorophyll.Concentration..OCI.Algorithm) %>%
  dplyr::select(-c(lat_dec, long_dec)) %>%
  mutate(extracted = "yes") # add column that specifies these ones are done
# this makes a big long list of chla estimates and the indiv paper number they're associated with 

# now make a dataframe of all the papers that don't have chla valuaes yet
meta_leftovers <- meta %>% 
  left_join(chla_prelim, by = "indiv") %>%
  replace_na(list(extracted = "no")) %>%
  filter(extracted == "no") %>% 
  dplyr::select(indiv, long_dec, lat_dec)

# Load chl a (yearly average for 2008)
raster_leftovers <- raster("Data/Rasters/A20021852021031.L3m_CU_CHL_chlor_a_4km.nc", warn = FALSE)

# turn the df into a shapefile
sf_leftovers <- meta_leftovers %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the chla
chla_leftovers <- raster::extract(raster_leftovers, sf_leftovers, df = TRUE) %>%
  bind_cols(meta_leftovers) %>%
  rename(chla_est = Chlorophyll.Concentration..OCI.Algorithm) %>% 
  dplyr::select(-c(lat_dec, long_dec, ID)) %>%
  mutate(extracted = "yes") #make it compatable to merge with earlier chla estimate df

# Add chla data to orig meta dataframe  -------
# now throw all the chla estimate data together, the specific years/seasons and the leftovers
chla_all_values <- rbind(chla_prelim, chla_leftovers)

meta_final <- meta %>%
  left_join(chla_all_values, by = "indiv")
# THERE ARE STILL SOME NA CHL A VALUES, MIGHT HAVE TO TWEAK THOSE COORDINATES

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
