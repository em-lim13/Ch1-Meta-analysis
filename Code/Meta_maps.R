# Making maps, baby
# July 27, 2023
# Em Lim

# Code to make a nice study site map in R

# Load packages ----
library(tidyverse)
library(sf)
library(maps)
library(ggplot2)
library(maptools) # wrld_simpl comes from here but is deprecated
library(rnaturalearth) # another simple world map?
library(ggspatial)

# Load meta_analysis data from Meta-analysis Script
meta_coords <- read_csv("Output/Data/lat_long2.csv") %>%
  as.data.frame() %>%
  mutate(lat = lat_dec,
         long = long_dec,
         climate = factor(climate, levels = c("Tropics", "Subtropics", "Temperate", "Polar"))) %>%
  st_as_sf(coords = c("long_dec", "lat_dec")) %>%
  st_set_crs(4326)

# Load simple world map from package "rnaturalearth"
world <- ne_countries(scale = "medium", returnclass = "sf")

# Simple map from maptools, but it's deprecated
data(wrld_simpl)
wrld_simpl <- wrld_simpl %>%
  st_as_sf() %>%
  st_transform(crs = 4326)

# OK SO
# When I plot the sp wrld_simpl using plot() the countries don't look squished
# But when I ggplot the sf version of wrld_simpl or the sf version of world (from ne_countries) that's when the countries look squished


# Good site for plotting help
# https://r-spatial.org/r/2018/10/25/ggplot2-sf.html


# Try plotting this ----
pal_ew <- c("#00496f", "#0f85a0", "#ed8b00", "#dd4124")
blue <- paste("#b9d1df", sep="")


ggplot() +
  geom_sf(data = world, fill = "white", colour = blue) +
  geom_sf(data = meta_coords, 
          aes(fill = climate), 
          colour = "black", pch = 21, size = 3, alpha = 0.5) +
  theme_bw() +
  theme(panel.background = element_rect(fill = blue)) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(expand = FALSE)+
  scale_y_continuous(limits = c(-79, 89)) +
  scale_x_continuous(limits = c(-180, 181)) +
  scale_fill_manual(values = (pal_ew))
  

# OLD shapefiles ----
# Load pieces of world map in 1:50 res
coast50 <- raster::shapefile("Shapefiles/50/ne_50m_coastline.shp") %>%
  st_as_sf() %>%
  st_transform(crs = 4326)

land50 <- raster::shapefile("Shapefiles/50/ne_50m_land.shp") %>%
  st_as_sf() %>%
  st_transform(crs = 4326)

ocean50 <- raster::shapefile("Shapefiles/50/ne_50m_ocean.shp") %>%
  st_as_sf() %>%
  st_transform(crs = 4326)

bound50 <- raster::shapefile("Shapefiles/50/ne_50m_admin_0_boundary_lines_land.shp") %>%
  st_as_sf() %>%
  st_transform(crs = 4326)

# Load pieces of world map in 1:110 res
coast110 <- raster::shapefile("Shapefiles/110/ne_110m_coastline.shp") %>%
  st_as_sf() %>%
  st_transform(crs = 4326)

land110 <- raster::shapefile("Shapefiles/110/ne_110m_land.shp") %>%
  st_as_sf() %>%
  st_transform(crs = 4326)

ocean110 <- raster::shapefile("Shapefiles/110/ne_110m_ocean.shp") %>%
  st_as_sf() %>%
  st_transform(crs = 4326)