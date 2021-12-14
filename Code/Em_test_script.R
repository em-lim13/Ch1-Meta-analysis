# Writing  a script to play with meta data and sea surface chl a data

# Seasonal chl a data
# https://oceancolor.gsfc.nasa.gov/cgi/l3?per=SN&prd=CHL_chlor_a.nc&sen=T&res=9km&num=24&ctg=Standard&date=21Dec1999

# Averaged chl a data 
# https://data.unep-wcmc.org/datasets/37



#install.packages('raster')
#install.packages('rgdal')
#install.packages('rasterVis')

#library(raster)
library(rgdal)
library(ggplot2)
library(ggmap)
library(rasterVis)
library(sf)
library(sp)
library(raster)
library(tidyverse)
library(PNWColors)

# Load in raster file.

global_chlra <- raster("Data/MeanSeaSurfaceChlorophyllA2015/01_Data/compositechla.ovr")

# load study sites
meta <- read_csv("Output/Data/metadata_pilot_chlra_extract.csv") %>%
  select(-c(Latitude, Longitude, continent, paper_num, obs_num) ) %>%
  mutate(Latitude = lat_dec,
         Longitude = long_dec) %>%
  select(-c(lat_dec, long_dec))

meta_spatial <- st_as_sf(meta, coords = c("Longitude", "Latitude"), crs = 4326)


# can I make a map
plot(global_chlra) #yikes that looks so bad the coordinates are fucked

# Basic summary stats:

chlra_stats_min <- cellStats(global_chlra, min) # Extracts minimum cell value
chlra_stats_median <- cellStats(global_chlra, median) # Extracts median cell value
chlra_stats_max <- cellStats(global_chlra, max) # Extracts maximum cell value
chlra_stats_average <- cellStats(global_chlra, average) # Extracts average cell value

##

# Quick histogram of the chlorophyll-a values - might be cool to compare the 
# distribution to the distribution of points in your extracted studies

chlra_hist <- hist(global_chlra, breaks = 30)

hist_data <- data.frame(counts = chlra_hist$counts, breaks = chlra_hist$mids)

ggplot(hist_data, aes(x = breaks, y = counts)) +
  geom_bar(stat = 'identity', alpha = 0.8) +
  xlab('Chlorophyll-a Concentration (mg/m^3)') +
  ylab('Number of 9 km^2 pixels assigned') +
  scale_fill_gradient(low = 'green', high = 'purple')

##

# Fun level plot of the global distribution of chlorophyll-a concentrations

levelplot(global_chlra, layers = 1, margin = list(FUN = 'median'), contour=TRUE)


# play with data

n <- length(unique(meta$chl_a))
bay_pal <- pnw_palette("Bay", n, type = "continuous")

ggplot(meta, aes(Longitude, Latitude, colour = chl_a)) +
  geom_point(size = 3) +
  scale_color_gradientn(colours = bay_pal)  +
  xlab("Longitude") +
  ylab("Latitude") + 
  labs(colour="Chl a") +
  geom_hline(yintercept =0, linetype = "dashed") +
  ylim(c(-75, 75)) +
  theme_classic() +
  theme(text = element_text(size = 18, colour = "black"))

ggsave("Output/Figures/lat_chla_fig.png", device = "png",
       height = 9, width = 16, dpi = 400)



#bin

arac <- read_csv("for_arcgis.csv")


arc <- arac %>%
  mutate(
    avg_chl = ifelse(chl_a < 5, "very low", ifelse(
      chl_a <10, "low", ifelse(
        chl_a <15, "medium", ifelse(
          chl_a <20, "high", ifelse(
            chl_a <25, "very high", "super high"
          )
        )
      )
    ))
  )

write_csv(arc, "data_for_map.csv")
