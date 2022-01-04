# Code to analyze mean and SE for Allgeier 2020 paper
# January 3, 2022


# Load packages
library(tidyverse)

# Load data
data <- read_csv("Data/Allgeier2020.csv") %>% 
  filter(keep_growth == 1) %>%
  mutate(ID = tagID)

# calculate mean, sd, and n for each species and each treatment
control <- data %>%
  filter(treat == "a") %>%
  group_by(sp)%>%
  summarise(mean = mean(g.area), sd = sd(g.area), n = n()) 

fish <- data %>%
  filter(treat == "b") %>%
  group_by(sp) %>%
  summarise(mean = mean(g.area), sd = sd(g.area), n = n()) 

# Load cell density data
cells <- read_csv("Data/Allgeier2020density.csv") %>%
  left_join(data, by = "ID")

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
