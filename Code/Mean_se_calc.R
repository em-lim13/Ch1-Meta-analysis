# Code to analyze mean and SE for Allgeier 2020, Shantz 2022, and Healing 2024 papers
# Em Lim
# January 3, 2022 updated Sept 2024


# Load packages
library(tidyverse)


# Allgeier 2020 ----
# Load data
allgeier <- read_csv("Data/Mean_SE_calc_data/Allgeier2020.csv") %>% 
  filter(keep_growth == 1) %>%
  mutate(ID = tagID)

# calculate mean, sd, and n for each species and each treatment
control <- allgeier %>%
  filter(treat == "a") %>%
  group_by(sp)%>%
  summarise(mean = mean(g.area), sd = sd(g.area), n = n()) 

fish <- allgeier %>%
  filter(treat == "b") %>%
  group_by(sp) %>%
  summarise(mean = mean(g.area), sd = sd(g.area), n = n()) 

# Load cell density data
cells <- read_csv("Data/Mean_SE_calc_data/Allgeier2020density.csv") %>%
  left_join(allgeier, by = "ID")

control_cells <- cells %>%
  filter(treat == "a") %>%
  drop_na(cell.cm) %>%
  group_by(sp)%>%
  summarise(mean = mean(cell.cm), sd = sd(cell.cm), n = n()) 

fish_cells <- cells %>%
  filter(treat == "b") %>%
  drop_na(cell.cm) %>%
  group_by(sp) %>%
  summarise(mean = mean(cell.cm), sd = sd(cell.cm), n = n()) 


# Shantz 2023 -----

# Load data
shantz <- read_csv("Data/Mean_SE_calc_data/Shantz2023.csv")

# calc mean +- SE
shantz_summary <- shantz %>%
  filter(temp.treatment == "Ambient" & day == 34) %>%
  group_by(fish.present) %>%
  summarise(zoo_avg = mean(avg.zoox.cm2), 
            zoo_sd = sd(avg.zoox.cm2), 
            prot_avg = mean(mg.prot.per.g), 
            prot_sd = sd(mg.prot.per.g), 
            growth_avg = mean(grow.mg.cm.d), 
            growth_sd = sd(grow.mg.cm.d), 
            n = n()) 

# one measure per day
shantz_summary_per_day <- shantz %>%
  filter(temp.treatment == "Ambient") %>%
  group_by(fish.present, day) %>%
  summarise(zoo_avg = mean(avg.zoox.cm2), 
            zoo_sd = sd(avg.zoox.cm2), 
            growth_avg = mean(grow.mg.cm.d), 
            growth_sd = sd(grow.mg.cm.d), 
            n = n()) 

# growth data?
shantz_grow_end <- shantz %>%
  filter(temp.treatment == "Ambient" & day == 34) %>%
  group_by(fish.present) %>%
  summarise(grow_end_avg = mean(grow.mg.cm.d), 
            grow_end_sd = sd(grow.mg.cm.d), 
            n = n()) 

# mean of both growth rates
shantz_grow_avg <- shantz %>%
  filter(temp.treatment == "Ambient") %>%
  group_by(fish.present) %>%
  summarise(grow_all_avg = mean(grow.mg.cm.d, na.rm = T), 
            grow_all_sd = sd(grow.mg.cm.d, na.rm = T)) 

shantz_grow <- shantz_grow_avg %>%
  left_join(shantz_grow_end, by = "fish.present") # VERY SIMILAR

# Healing 2024 ------
healing <- read_csv("Data/mean_SE_calc_data/Healing2024.csv") 

healing_0 <- healing %>%
  select(-GPS) %>%
  replace(is.na(.), 0)

healing_long <- healing_0 %>%
  filter(Marker == "Average") %>%
  pivot_longer(cols = `Bladder_wrack_%`:Common_shore_crab,
               names_to = "species",
               values_to = "percent")

healing_sum <- healing_long %>%
  group_by(species, `Control/treatment`) %>%
  summarise(mean = mean(percent),
            sd = sd(percent),
            n = n())


