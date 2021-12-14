##

# Nutrient Cycling Lit Review Script
# Date: Aug 13 2019
# Author: Em Lim
# Email: em_lim@sfu.ca

# Description: Script using the package 
#'devtools' to manage and process exported
# .bib files retrieved from literature search.

# 1.) Install package as needed, load library:

#install.packages('revtools')
#install.packages('tidyverse')
#install.packages('devtools')
#install.packages('wosr')

library(devtools)
library(revtools)
library(tidyverse)
library(wosr)


# 2.) Load exported .bib/.ris file from database.

# Search terms: (ammon* OR nitrog* OR phosph* OR nutrient) 
# AND (excret* OR egest* OR mutual* OR symbio* OR subsid* 
# OR consumer-mediat*) AND (coral* OR kelp* OR alga* 
# OR macroalgae OR microalgae OR phytoplankton OR zooxanthell*
# OR primary prod* OR seagrass* OR sea grass* OR seaweed* 
# OR mangrove*) AND (reef OR subtidal OR intertidal OR coast* 
# OR marine OR benthic OR coral OR kelp) AND (grow* OR grew 
# OR biomass OR abundance OR weight)

# Web of Science search, August 13, 2020
# TS= ((ammon* OR nitrog* OR phosph* OR nutrient) AND (excret* OR egest* OR mutual* OR symbio* OR subsid* OR consumer-mediat*) AND (coral* OR kelp* OR alga* OR macroalgae OR microalgae OR phytoplankton OR zooxanthell* OR primary prod* OR seagrass* OR sea grass* OR seaweed* OR mangrove*) AND (reef OR subtidal OR intertidal OR coast* OR marine OR benthic OR coral OR kelp) AND (grow* OR grew OR biomass OR abundance OR weight))
# 1,215 hits

# Refined by: DOCUMENT TYPES: ( ARTICLE OR PROCEEDINGS PAPER )
# 1,167

# Refined by: WEB OF SCIENCE CATEGORIES: ( MARINE FRESHWATER BIOLOGY OR ENGINEERING ENVIRONMENTAL OR ECOLOGY OR OCEANOGRAPHY OR ZOOLOGY OR ENVIRONMENTAL SCIENCES OR LIMNOLOGY OR MULTIDISCIPLINARY SCIENCES OR BIOLOGY OR PLANT SCIENCES OR FISHERIES OR BIODIVERSITY CONSERVATION )
# 1006 hits

setwd("~/Documents/MSc/Directed Readings/Meta_Analysis/Bib files")

WoSa <- read_bibliography('wos_search_Aug2020a.bib')
WoSb <- read_bibliography('wos_search_Aug2020b.bib')
WoSc <- read_bibliography('wos_search_Aug2020c.bib')

WoS <- bind_rows(WoSa, WoSb, WoSc)

# ASFA search, August 13 2020: 

# noft((ammon* OR nitrog* OR phosph* OR nutrient) AND (excret* OR egest* OR mutual* OR symbio* OR subsid* OR consumer-mediat*) AND (coral* OR kelp* OR alga* OR macroalgae OR microalgae OR phytoplankton OR zooxanthell* OR primary prod* OR seagrass* OR sea grass* OR seaweed* OR mangrove*) AND (reef OR subtidal OR intertidal OR coast* OR marine OR benthic OR coral OR kelp) AND (grow* OR grew OR biomass OR abundance OR weight)) 
# Hits: 1428

# Limit to Scholarly type, reports
# Hits: 1167

# NOT mainsubject(freshwater OR bacteria OR aquaculture OR lakes OR rivers OR cyanobacteria OR biogeochemistry OR nitrogen fixation)
# Hits: 784

ASFA <- read_bibliography('asfa_search_Aug2020.ris')


# combine both searches
combo_WoS_ASFA <- bind_rows(WoS, ASFA)
# Hits: 1790

# 3.) Identify and remove duplicates based on title:

combo_match <- find_duplicates(combo_WoS_ASFA, match_variable = "title", group_variables = NULL, match_function = "fuzzdist",
                             method = "fuzz_m_ratio", threshold = 0.1, to_lower = TRUE, remove_punctuation = FALSE)

combo_unique <- extract_unique_references(combo_WoS_ASFA, combo_match)

save(combo_unique, file = "combo_unique.rData")

#load("combo_unique.rData")

# 1429 entries once duplicates removed

# 4.) Screen titles and abstracts using 'revtools' GUI:

# Criteria for inclusion/exclusion:

# [1.] Marine environment                             Y/N, yes included
# [2.] Animals are the source of nutrients?           Y/N, yes included
# [3.] Primary producer measured in response          Y/N, yes included


combo_title_screened <- screen_titles(combo_unique) 


#make the screened file into a bib
#write_bibliography(combo_title_screened, 'combo_titles_screened_2020.bib', format = 'bib')

# Can jump forward to here if you need to load data -------
#import that bib
all_title_screened <- read_bibliography('Bib_files/combo_titles_screened_2020.bib')

#remove the excluded titles
combo_subset <- all_title_screened[all_title_screened$screened_titles != 'excluded',]

# Screening the titles --> 199 hits

# Screen the abstracts for inclusion, same criteria

#abstracts_screened <- screen_abstracts(combo_subset)

#write_bibliography(abstracts_screened, 'abstracts_screened_2020.bib', format = 'bib')
abstracts_screened <- read_bibliography("abstracts_screened_2020.bib")

abstracts_screened_subset <- abstracts_screened[abstracts_screened$screened_abstracts != 'excluded',]

#write.csv(abstracts_screened_subset, 'final_abstracts.csv')

titles <- abstracts_screened_subset %>% select(title)

# screening abstracts = 100 hits

# 5.) Concatenate final list of papers from ASFA and WoS, filter for duplicates and retrieve 
# records for final review list.


write.csv(titles, 'final_titles.csv')

## papers_final <- merge(x = titles_unique, y = ASFA_abstract_screened_subset, by = 'title', all.y = TRUE)
