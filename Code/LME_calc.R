# Large Marine Ecoregions
# Em Lim
# try to get ecoregions

# load packages
library(tidyverse)
library(ggspatial)
library(sf)
sf_use_s2(FALSE)


# load data
lme <- sf::st_read("Data/lme_shp/lme66.shp") %>%
  st_sf() %>%
  st_set_crs(4326)

# offshore <- sf::st_read("Data/offshore_shp/lme66-offshore.shp") %>%
#   st_sf() %>%
#   st_set_crs(4326)

# load coords
coords <- read_csv("Data/Output_data/meta_df.csv") %>%
  distinct(paper_num, title, lat_dec, long_dec) %>%
  as.data.frame() %>%
  mutate(lat = case_when(
    paper_num == "24" ~ -45.734127, 
    paper_num == "25" ~ 44.275933, 
    paper_num == "30" ~ 48.261035, 
    paper_num == "32" ~ 17.839382, 
    paper_num == "33" ~ 48.438672, 
    paper_num == "35" ~ 26.255705, 
    paper_num == "37" ~ 17.606429, 
    paper_num == "38" ~ 33, 
    paper_num == "43" ~ 40.585689, 
    paper_num == "45" ~ 32.729070, 
    
    TRUE ~ as.numeric(lat_dec)),
         long = 
      case_when(
        paper_num == "24" ~ 170.869691, 
        paper_num == "25" ~ -124.287282, 
        paper_num == "30" ~ -4.945408, 
        paper_num == "32" ~ -64.594842, 
        paper_num == "33" ~ -124.80177, 
        paper_num == "35" ~ -77.6, 
        paper_num == "37" ~ -64.675044, 
        paper_num == "38" ~ 32, 
        paper_num == "43" ~ 14.017679, 
        paper_num == "45" ~ -117.312799, 
        
        TRUE ~ as.numeric(long_dec))
    ) %>%
  st_as_sf(coords = c("long", "lat")) %>%
  st_set_crs(4326)


# try to join up points
coords_LME <- st_join(coords, lme %>% select(LME_NUMBER, LME_NAME)) %>%
  mutate(LME_NAME = case_when(
    paper_num %in% c("1", "17","21", "39", "41", "42") ~ "Offshore",
    TRUE ~ as.character(LME_NAME)
  ))

# missing points
offshore <- meta_df %>%
  filter(paper_num %in% c("1", "17","21", "39", "41", "42")) %>%
  select(paper_num, country, location, continent, ocean, region, lat_dec, long_dec)%>%
  unique() %>%
  st_as_sf(coords = c("long_dec", "lat_dec")) %>%
  st_set_crs(4326)

# try to plot these missing offshore points
ggplot() +
  geom_sf(data = lme, aes(fill = LME_NAME), color = "black", alpha = 0.6) +
  geom_sf(data = offshore, color = "red", size = 2) +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(legend.position = "none")

# save the LME classifications
lme_coords <- coords_LME %>%
  select(paper_num, LME_NAME, LME_NUMBER) %>%
  st_drop_geometry() %>%
  unique() %>%
  rename(lme_name = LME_NAME,
         lme_num = LME_NUMBER)
# write_csv(lme_coords, "Data/Output_data/LME_df.csv")

lme_list <- lme %>%
  st_drop_geometry() %>%
  select(LME_NUMBER, LME_NAME) %>%
  arrange(desc(-LME_NUMBER)) %>%
  unite(full_name, LME_NUMBER, LME_NAME, sep = ". ")

#write_csv(lme_list, "Tables/LME_list.csv")  
