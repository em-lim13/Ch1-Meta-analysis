# SST again

# laod packages
library(dplyr) # A staple for modern data management in R
library(tidyverse)
library(lubridate) # Useful functions for dealing with dates
library(ggplot2) # The preferred library for data visualisation
library(tidync) # For easily dealing with NetCDF data
library(rerddap) # For easily downloading subsets of data
library(doParallel) # For parallel processing

# help website here:
# https://cran.r-project.org/web/packages/heatwaveR/vignettes/OISST_preparation.html

# The information for the NOAA OISST data
rerddap::info(datasetid = "ncdcOisst21Agg_LonPM180", url = "https://coastwatch.pfeg.noaa.gov/erddap/")

# Note that there is also a version with lon values from 0 yo 360
rerddap::info(datasetid = "ncdcOisst21Agg", url = "https://coastwatch.pfeg.noaa.gov/erddap/")


# Load meta-analysis data
meta <- read_csv("Data/Output_data/meta_df.csv") %>%
  # make both year and year of study character
  mutate(year = as.character(year),
         year_of_study = if_else(is.na(year_of_study), year, year_of_study)) %>%
  # just select and rename needed columns
  select(indiv, lat_dec, long_dec, season, year_of_study) %>%
  rename(year = year_of_study) %>%
  # spread year ranges so each year has it's own row
  separate(year, into = c("start", "end"), sep = "-", fill = "right") %>%
  mutate(across(c(start, end), as.integer),
         end = if_else(is.na(end), start, end)) %>%
  mutate(
    start_date = as.Date(paste0(start, "-01-01")),
    end_date   = as.Date(paste0(end, "-12-31")),
    # shift start and end dates for periods before the sst data started
    start_date = case_when(
      start_date <= as.Date("1981-01-01") ~ as.Date("1982-01-01"),
      TRUE ~ start_date),
    end_date = case_when(
      end_date <= as.Date("1981-12-31") ~ as.Date("1982-12-31"),
      TRUE ~ end_date))

# make df of unique requests
unique_requests <- meta %>%
  distinct(lat_dec, long_dec, start_date, end_date) %>%
  mutate(
    request_id = row_number(),
    lat_dec  = ifelse(request_id == 64, -40.937904, lat_dec),
    long_dec = ifelse(request_id == 64, 174.137357, long_dec)
  )

# download SST data using original function (this will take 30 - 60 mins)
# unhash to run, hashed out so i don't run by accident
# sst_downloaded <- unique_requests %>% # function written in fuction script
#   split(seq_len(nrow(.))) %>%
#   map_dfr(download_sst)

q <- sst_downloaded %>%
  count(request_id)

# missing 2, 41, 42, 52, 53, 57, 58, 63, 66
# 52, 53, 57, 58 out of range
# I think the rest aren't in a cell with a value! coordinates in inlet or on land

# Run fancy download function for the NA cells -------

# filter unique dataframe down to just the missed cells
unique_missing <- unique_requests %>%
  filter(!(request_id %in% sst_downloaded$request_id)) 

# download data for those originally NA cells
sst_downloaded_m <- unique_missing %>%
  split(seq_len(nrow(.))) %>%
  map_dfr(download_sst_nearest)

# Put all of the data together
sst_downloaded_full <- rbind(sst_downloaded, sst_downloaded_m)

# clean the df to save it
sst_downloaded_clean <- sst_downloaded_full %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    lat_dec = as.numeric(lat_dec),
    long_dec = as.numeric(long_dec),
    t = as.Date(t),   
    temp = as.numeric(temp)
  ) %>%
  select(request_id, lat_dec, long_dec, start_date, end_date, lon, lat, t, temp)

# write_csv(sst_downloaded_full, "Data/Output_data/sst_data_last_one.csv")

# Jump here to work with already downloaded SST data ----------
sst_downloaded_clean <- read_csv("Data/Output_data/sst_data_last_one.csv")

# Get annual averages
annual_sst <- sst_downloaded_clean %>%
  group_by(request_id, lat_dec, long_dec, start_date, end_date) %>%
  summarise(annual_sst = mean(temp, na.rm = TRUE)) %>%
  mutate(
    lat_dec  = ifelse(request_id == 64, -41.151667, lat_dec),
    long_dec = ifelse(request_id == 64, 173.875000, long_dec)
  )


# join the sst data back to the meta df!
meta_sst <- meta %>%
  left_join(annual_sst, by = c("lat_dec", "long_dec", "start_date", "end_date")) %>%
  select(indiv, annual_sst)

# write_csv(meta_sst, "Data/Output_data/meta_sst_final.csv")