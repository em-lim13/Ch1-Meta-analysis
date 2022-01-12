# Script to extract Chl a data for the meta analysis
# Em Lim
# Jan 4, 2021

# LOOK AT DIFFERENT SENSORS BC THEY COVER DIFFERENT TIME PERIODS!!!!!!

# Be warned what follows is a bit jank
# Be careful with this..... 
# The raster files are currently living on my harddrive because they are too big for my laptop or github
# If you need to rerun this, move that raster file into the data folder
# But the word doc still has all the download links

# Load packages and data -------
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
library(tidyverse)
library(metafor)

# Load data (THIS IS FROM META-ANALYSIS SCRIPT)
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

# Chl a data extraction------

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


# WINTER 1996 -----
# Load a chl a file
# chl a value is in mg/m^2

raster_winter1996 <- raster("Data/Rasters/O19963561997079.L3m_SNWI_CHL_chlor_a_9km.nc", warn = FALSE)

# make a df for this one season / year
meta_winter1996 <- meta %>%
  filter(season == "winter" & year_of_study == "1996") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_winter1996 <- meta_winter1996 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_winter1996 <- raster::extract(raster_winter1996, sf_winter1996, df = TRUE) %>%
  bind_cols(meta_winter1996) %>%
  dplyr::select(-ID)

# SPRING 1997 -----
# Load a chl a file
raster_spring1997 <- raster("Data/Rasters/O19970801997171.L3m_SNSP_CHL_chlor_a_9km.nc", warn = FALSE)

# make a df for this one season / year
meta_spring1997 <- meta %>%
  filter(season == "spring" & year_of_study == "1997") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_spring1997 <- meta_spring1997 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_spring1997 <- raster::extract(raster_spring1997, sf_spring1997, df = TRUE) %>%
  bind_cols(meta_spring1997) %>%
  dplyr::select(-ID)

# SUMMER 1997 -----
# Load a chl a file
raster_summer1997 <- raster("Data/Rasters/S19971721997263.L3m_SNSU_CHL_chlor_a_9km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer1997 <- meta %>%
  filter(season == "summer" & year_of_study == "1997") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer1997 <- meta_summer1997 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer1997 <- raster::extract(raster_summer1997, sf_summer1997, df = TRUE) %>%
  bind_cols(meta_summer1997) %>%
  dplyr::select(-ID)

# FALL 1997 -----
# Load a chl a file
raster_fall1997 <- raster("Data/Rasters/S19972641997354.L3m_SNAU_CHL_chlor_a_9km.nc", warn = FALSE)

# make a df for this one season / year
meta_fall1997 <- meta %>%
  filter(season == "fall" & year_of_study == "1997") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_fall1997 <- meta_fall1997 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_fall1997 <- raster::extract(raster_fall1997, sf_fall1997, df = TRUE) %>%
  bind_cols(meta_fall1997) %>%
  dplyr::select(-ID)

# SPRING 2002 -----
# Load a chl a file
# chl a value is in mg/m^2
# using the summer 2002 file bc there isn't a spring one and it's close enough
raster_spring2002 <- raster("Data/Rasters/A20021722002263.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_spring2002 <- meta %>%
  filter(season == "spring" & year_of_study == "2002") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_spring2002 <- meta_spring2002 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_spring2002 <- raster::extract(raster_spring2002, sf_spring2002, df = TRUE) %>%
  bind_cols(meta_spring2002) %>%
  dplyr::select(-ID)


# SUMMER 2004 -----
# Load a chl a file
raster_summer2004 <- raster("Data/Rasters/A20041732004264.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2004 <- meta %>%
  filter(season == "summer" & year_of_study == "2004") %>%
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
# Try using summer 2016 raster since summer 2015 isn't working
raster_summer2015 <- raster("Data/Rasters/A20161732016264.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

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

# SPRING 2018 -------
# Load chl a
raster_spring2018 <- raster("Data/Rasters/A20180802018171.L3m_SNSP_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_spring2018 <- meta %>%
  filter(season == "spring" & year_of_study == "2018") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_spring2018 <- meta_spring2018 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_spring2018 <- raster::extract(raster_spring2018, sf_spring2018, df = TRUE) %>%
  bind_cols(meta_spring2018) %>%
  dplyr::select(-ID)

# SUMMER 2018 -------
# Load chl a
raster_summer2018 <- raster("Data/Rasters/A20181722018263.L3m_SNSU_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_summer2018 <- meta %>%
  filter(season == "summer" & year_of_study == "2018") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_summer2018 <- meta_summer2018 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the 
chla_summer2018 <- raster::extract(raster_summer2018, sf_summer2018, df = TRUE) %>%
  bind_cols(meta_summer2018) %>%
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

# 1979 - 1981 -------
# Load chl a from 1979, 1980, 1981
raster_1979 <- raster("Data/Rasters/C19790011979365.L3m_YR_CHL_chl_ocx_4km.nc", warn = FALSE)

raster_1980 <- raster("Data/Rasters/C19800011980366.L3m_YR_CHL_chl_ocx_4km.nc", warn = FALSE)

raster_1981 <- raster("Data/Rasters/C19810011981365.L3m_YR_CHL_chl_ocx_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_1979_1981 <- meta %>%
  filter(year_of_study == "1979-1981") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_1979_1981 <- meta_1979_1981 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the chl a data
temp_1979 <- raster::extract(raster_1979, sf_1979_1981, df = TRUE) %>%
  rename(chla_1979 = Chlorophyll.Concentration..OC3.Algorithm )%>%
  dplyr::select(-ID)

temp_1980 <- raster::extract(raster_1980, sf_1979_1981, df = TRUE) %>%
  rename(chla_1980 = Chlorophyll.Concentration..OC3.Algorithm )%>%
  dplyr::select(-ID)

temp_1981 <- raster::extract(raster_1981, sf_1979_1981, df = TRUE) %>%
  rename(chla_1981 = Chlorophyll.Concentration..OC3.Algorithm )%>%
  dplyr::select(-ID)

# now bind all these yearly averages into the same dataframe and average across all years
chla_1979_1981 <- bind_cols(temp_1979, temp_1980, temp_1981) %>%
  mutate(Chlorophyll.Concentration..OCI.Algorithm = rowMeans(cbind(chla_1979, chla_1980, chla_1981))) %>%
  dplyr::select(Chlorophyll.Concentration..OCI.Algorithm) %>%
  bind_cols(meta_1979_1981)
# the chl a name is different but that will make it compatible with the rest

# 1998 - 2000 --------
# Load chl a from 1998, 1999, 2000
raster_1998 <- raster("Data/Rasters/S19980011998365.L3m_YR_CHL_chlor_a_9km.nc", warn = FALSE)

raster_1999 <- raster("Data/Rasters/S19990011999365.L3m_YR_CHL_chlor_a_9km.nc", warn = FALSE)

raster_2000 <- raster("Data/Rasters/S20000012000366.L3m_YR_CHL_chlor_a_9km.nc", warn = FALSE)

# make a df for this one season / year
meta_1998_2000 <- meta %>%
  filter(year_of_study == "1998-2000") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_1998_2000 <- meta_1998_2000 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the chl a data
temp_1998 <- raster::extract(raster_1998, sf_1998_2000, df = TRUE) %>%
  rename(chla_1998 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_1999 <- raster::extract(raster_1999, sf_1998_2000, df = TRUE) %>%
  rename(chla_1999 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_2000 <- raster::extract(raster_2000, sf_1998_2000, df = TRUE) %>%
  rename(chla_2000 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

# now bind all these yearly averages into the same dataframe and average across all years
chla_1998_2000 <- bind_cols(temp_1998, temp_1999, temp_2000) %>%
  mutate(Chlorophyll.Concentration..OCI.Algorithm = rowMeans(cbind(chla_1998, chla_1999, chla_2000))) %>%
  dplyr::select(Chlorophyll.Concentration..OCI.Algorithm) %>%
  bind_cols(meta_1998_2000)

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

# 2005 - 2006 -------
# Load chl a from 2005, 2006
raster_2005 <- raster("Data/Rasters/A20050012005365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

raster_2006 <- raster("Data/Rasters/A20060012006365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_2005_2006 <- meta %>%
  filter(year_of_study == "2005-2006") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_2005_2006 <- meta_2005_2006 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the chl a data
temp_2005b <- raster::extract(raster_2005, sf_2005_2006, df = TRUE) %>%
  rename(chla_2005 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_2006 <- raster::extract(raster_2006, sf_2005_2006, df = TRUE) %>%
  rename(chla_2006 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

# now bind all these yearly averages into the same dataframe and average across all years
chla_2005_2006 <- bind_cols(temp_2005b, temp_2006) %>%
  mutate(Chlorophyll.Concentration..OCI.Algorithm = rowMeans(cbind(chla_2005, chla_2006))) %>%
  dplyr::select(Chlorophyll.Concentration..OCI.Algorithm) %>%
  bind_cols(meta_2005_2006)

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

# 2011 - 2013 -------
# Load yearly average chl a from 2013, 2014
raster_2011 <- raster("Data/Rasters/A20110012011365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

raster_2012 <- raster("Data/Rasters/A20120012012366.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

raster_2013 <- raster("Data/Rasters/A20130012013365.L3m_YR_CHL_chlor_a_4km.nc", warn = FALSE)

# make a df for this one season / year
meta_2011_2013 <- meta %>%
  filter(year_of_study == "2011-2013") %>%
  dplyr::select(indiv, long_dec, lat_dec)

# turn the df into a shapefile
sf_2011_2013 <- meta_2011_2013 %>%
  dplyr::select(-indiv) %>%
  st_as_sf(coords = c("long_dec", "lat_dec"), crs = my.projection) %>%
  as(Class = "Spatial")

# use the raster and the shapefile to extract the chl a data
temp_2011 <- raster::extract(raster_2011, sf_2011_2013, df = TRUE) %>%
  rename(chla_2011 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_2012 <- raster::extract(raster_2012, sf_2011_2013, df = TRUE) %>%
  rename(chla_2012 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

temp_2013b <- raster::extract(raster_2013, sf_2011_2013, df = TRUE) %>%
  rename(chla_2013 = Chlorophyll.Concentration..OCI.Algorithm )%>%
  dplyr::select(-ID)

# now bind all these yearly averages into the same dataframe and average across all years
chla_2011_2013 <- bind_cols(temp_2011, temp_2012, temp_2013b) %>%
  mutate(Chlorophyll.Concentration..OCI.Algorithm = rowMeans(cbind(chla_2011, chla_2012, chla_2013))) %>%
  dplyr::select(Chlorophyll.Concentration..OCI.Algorithm) %>%
  bind_cols(meta_2011_2013)

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

# Combine chl_a dataframes -------

# Stick all the completed chla dataframes together
chla_prelim <- rbind(chla_1979_1981, chla_1998_2000, chla_2000_2005, chla_2005_2006, chla_2002, chla_2003_2005, chla_2008, chla_2011_2013, chla_2013_2014, chla_2015_2016, chla_fall1997, chla_fall2012, chla_fall2019, chla_spring1997, chla_spring2002, chla_spring2005, chla_spring2015, chla_spring2016, chla_spring2018, chla_summer1997, chla_summer2004, chla_summer2005, chla_summer2006, chla_summer2009, chla_summer2010, chla_summer2011, chla_summer2012, chla_summer2013, chla_summer2015, chla_summer2018, chla_winter1996, chla_winter2005, chla_winter2010, chla_winter2013, chla_winter2017) %>% #24 dfs 
  rename(chla_est = Chlorophyll.Concentration..OCI.Algorithm) %>%
  dplyr::select(-c(lat_dec, long_dec)) %>%
  mutate(extracted = "yes") # add column that specifies these ones are done
# this makes a big long list of chla estimates and the indiv paper number they're associated with 

# LEFTOVERS ------
# maybe try to find a better match for the misc papers?
# Eg if it's a 1995 paper, could use a summer rolling avg

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


# now throw all the chla estimate data together, the specific years/seasons and the leftovers
chla_all_values <- rbind(chla_prelim, chla_leftovers)
# will be "missing some" bc I didn't look at the papers with negative effect sizes

#write_csv(chla_all_values, path =  "Output/Data/all_extracted_chla.csv")
# write a csv that I can pick up in the Meta-analysis Script


# Papers with NA's: 

# Summer 1997


# ALL FIXED NOW

# Hurd et al 2000
  # summer 1995 -> in the leftover extraction
  # Neither Dixon or Kiixin worked
  # The whole Barkley sound region has INCREDIBLY high chl a 
  # set both closer to middle of trevor channel

# Chase et al 2018:
# The spring 2016 coordinates worked fine, but the summer 2015 coordinates didn't work
# Not sure if that season/year is just spotty with coverage?
# See if I change summer 2015 to summer 2016?????

# Bracken 2004
  # spring 2002
  # Inside a bay, successfully solved the problem by moving the coordinates just outside the bay

# Archer et al 2015
# Summer 2013
# Sort of in a bay too
# Tried moving the coordinates to the right and that worked
