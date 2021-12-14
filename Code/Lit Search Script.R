# Nutrient Cycling Lit Review Script
# Date: Aug 13 2019
# Author: Em Lim

# Description: Script using the package 'devtools' to manage and process exported .bib files retrieved from literature search.

# Load packages -------
library(devtools)
library(revtools)
library(tidyverse)
library(wosr)
library(synthesisr)

# 1900 - 2020/08/13 -------

# 1) Load exported .bib/.ris file from database ----

# Web of Science search, August 13, 2020
# TS= ((ammon* OR nitrog* OR phosph* OR nutrient) AND (excret* OR egest* OR mutual* OR symbio* OR subsid* OR consumer-mediat*) AND (coral* OR kelp* OR alga* OR macroalgae OR microalgae OR phytoplankton OR zooxanthell* OR primary prod* OR seagrass* OR sea grass* OR seaweed* OR mangrove*) AND (reef OR subtidal OR intertidal OR coast* OR marine OR benthic OR coral OR kelp) AND (grow* OR grew OR biomass OR abundance OR weight))
# 1,215 hits (redone in 2024 = 1,325 hits)

# Refined by: DOCUMENT TYPES: ( ARTICLE OR PROCEEDINGS PAPER )
# 1,167 (redone in 2024 = 1,222 hits)

# Refined by: WEB OF SCIENCE CATEGORIES: ( MARINE FRESHWATER BIOLOGY OR ENGINEERING ENVIRONMENTAL OR ECOLOGY OR OCEANOGRAPHY OR ZOOLOGY OR ENVIRONMENTAL SCIENCES OR LIMNOLOGY OR MULTIDISCIPLINARY SCIENCES OR BIOLOGY OR PLANT SCIENCES OR FISHERIES OR BIODIVERSITY CONSERVATION )
# 1006 hits (redone in 2024 = 1067 hits)

WoS <- bind_rows(
  read_bibliography('Data/Bib_files/Input/wos_search_Aug2020a.bib'),
  read_bibliography('Data/Bib_files/Input/wos_search_Aug2020b.bib'),
  read_bibliography('Data/Bib_files/Input/wos_search_Aug2020c.bib'))

# ASFA search, August 13 2020: 
# Aquatic Sciences and Fisheries Abstracts

# noft((ammon* OR nitrog* OR phosph* OR nutrient) AND (excret* OR egest* OR mutual* OR symbio* OR subsid* OR consumer-mediat*) AND (coral* OR kelp* OR alga* OR macroalgae OR microalgae OR phytoplankton OR zooxanthell* OR primary prod* OR seagrass* OR sea grass* OR seaweed* OR mangrove*) AND (reef OR subtidal OR intertidal OR coast* OR marine OR benthic OR coral OR kelp) AND (grow* OR grew OR biomass OR abundance OR weight)) 
# Hits: 1428 (redone in 2024 = 1,432)

# Limit to Scholarly type, reports
# Hits: 1167 (redone in 2024 = 1,177)

# NOT mainsubject(freshwater OR bacteria OR aquaculture OR lakes OR rivers OR cyanobacteria OR biogeochemistry OR nitrogen fixation)
# Hits: 784 (redone in 2024 = 828)

ASFA <- synthesisr::read_ref('Data/Bib_files/Input/asfa_search_Aug2020.ris')

# combine WoS and ASFA searches 
combo_WoS_ASFA <- bind_rows(WoS, ASFA)
# Hits: 1790

# 2) Identify and remove duplicates based on title: ------
combo_match <- revtools::find_duplicates(combo_WoS_ASFA, match_variable = "title", group_variables = NULL, match_function = "fuzzdist", method = "fuzz_m_ratio", threshold = 0.1, to_lower = TRUE, remove_punctuation = FALSE)

combo_unique <- extract_unique_references(combo_WoS_ASFA, combo_match)
# 1429 entries once duplicates removed

# Save this df so I can just jump back in here
#save(combo_unique, file = "Data/Bib_files/Output/combo_unique.rData")
#load("Data/Bib_filesOutput/combo_unique.rData")


# 3) Screen titles and abstracts using 'revtools' GUI: --------

# Criteria for inclusion/exclusion:
# [1] Marine environment                             Y/N, yes included
# [2] Animals are the source of nutrients?           Y/N, yes included
# [3] Primary producer measured in response          Y/N, yes included

combo_title_screened <- screen_titles(combo_unique) # this opens the app to screen titles

#make the screened file into a bib
#write_bibliography(combo_title_screened, 'Data/Bib_files/Output/combo_titles_screened_2020.bib', format = 'bib')

#import that bib
all_title_screened <- read_bibliography('Data/Bib_files/Output/combo_titles_screened_2020.bib')

#remove the excluded titles
combo_subset <- all_title_screened %>%
  filter(screened_titles != "excluded")
# Screening the titles --> 199 hits

# Screen the abstracts for inclusion, same criteria
abstracts_screened <- screen_abstracts(combo_subset) # this opens app

# write the screened abstracts to a file
#write_bibliography(abstracts_screened, 'Data/Bib_files/Output/abstracts_screened_2020.bib', format = 'bib')

abstracts_screened <- read_bibliography("Data/Bib_files/Output/abstracts_screened_2020.bib")

abstracts_screened_subset <- abstracts_screened %>%
  filter(screened_abstracts != "excluded")

#write.csv(abstracts_screened_subset, 'Data/Bib_files/Output/final_abstracts.csv')

titles <- abstracts_screened_subset %>% select(title)

# screening abstracts = 100 hits

# 5) Concatenate final list of papers from ASFA and WoS, filter for duplicates and retrieve records for final review list.

# write.csv(titles, 'Data/Bib_files/Output/final_titles.csv')



# Redo the 1900 - 2020 search -----

# reload the 1900 - 2020 wos search I did in 2024
wos_redo <- bind_rows(synthesisr::read_ref('Data/Bib_files/Input/2024/wos_search_Aug2020_redone_a.ris'), synthesisr::read_ref('Data/Bib_files/Input/2024/wos_search_Aug2020_redone_b.ris'))

# load redone 1900 - 2020 ASFA search
asfa_redo <- synthesisr::read_ref('Data/Bib_files/Input/2024/asfa_search_Aug2020_redone.ris')

# combine redone wos and asfa searches
combo_wos_asfa_redo <- bind_rows(wos_redo, asfa_redo)

# 2) Identify and remove duplicates based on title:
combo_match_redo <- revtools::find_duplicates(combo_wos_asfa_redo, match_variable = "title", group_variables = NULL, match_function = "fuzzdist", method = "fuzz_m_ratio", threshold = 0.1, to_lower = TRUE, remove_punctuation = FALSE)

# now just list those unique hits
combo_unique_redo <- revtools::extract_unique_references(combo_wos_asfa_redo, combo_match_redo)
# 1571 unique refs, 324 duplicates

# OK and NOW let's see how many "new" hits redoing the search got

# take out italics characters bc those seem to mess things up
combo_unique_redo_fixed <- combo_unique_redo %>%
  mutate(title = gsub('<i>', '', title),
         title = gsub('</i>', '', title))

combo_unique_fixed <- combo_unique %>%
  mutate(title = gsub('<i>', '', title),
         title = gsub('</i>', '', title))

redo_2020 <- combo_unique_redo_fixed %>%
  filter(!(title %in% combo_unique_fixed$title)) # 178 hits

# this opens the app to screen titles
# title_screened_redo <- screen_titles(redo_2020) 

#make the screened file into a bib
# write_bibliography(title_screened_redo, 'Data/Bib_files/Output/title_screened_redo.bib', format = 'bib')

#import that bib
all_title_screened_redo <- synthesisr::read_ref('Data/Bib_files/Output/title_screened_redo.bib')

#remove the excluded titles
combo_subset_redo <- all_title_screened_redo %>%
  filter(screened_titles != "excluded") %>%
  mutate(label = 1:n())

# Screening the titles --> 24 hits
# some of these were definitely in the original search

# Screen the abstracts for inclusion, same criteria
# abstracts_screened_redo <- screen_abstracts(combo_subset_redo) # this opens app

# write the screened abstracts to a file
# write_bibliography(abstracts_screened_redo, 'Data/Bib_files/Output/abstracts_screened_redo_2020.bib', format = 'bib')

abstracts_screened_redo <- revtools::read_bibliography("Data/Bib_files/Output/abstracts_screened_redo_2020.bib")

abstracts_screened_subset_redo <- abstracts_screened_redo %>%
  filter(screened_abstracts != "excluded") %>% # 9 hits, double check these aren't already in the sheet
  filter(label != "X15" ) %>%
  filter(label != "X21") %>%
  filter(label != "X6") %>%
  filter(label != "X7") %>%
  filter(label != "X10") %>%
  filter(label != "X24") %>%
  filter(label != "X17") # 2 once papers that were in the original search were removed
  # I don't think I can actually include either of those papers due to a lack of controls/direct measurement of pp as function of cnc

# write.csv(abstracts_screened_subset_redo, 'Data/Bib_files/Output/final_titles_redo.csv')


# 2020/08/14 - 2022/01/03 ------

# 1) Load .ris files

# Web of Science search, August 13 2020 - January 3, 2022

# Same search string as in 2020
# 36 hits (redone in 2024 = 139 hits)

# Refined by: DOCUMENT TYPES: ( ARTICLE OR PROCEEDINGS PAPER )
# 35 hits (Redone in 2024 = 232 hits)

# Refined by same categories as in 2020
# 26 hits (Redone in 2024 = 101 hits)

WoS2022 <- synthesisr::read_ref('Data/Bib_files/Input/wos_search_Jan2022.ris')
# this is the new WoS list


# AGAIN for ASFA
# ASFA search, January 3 2022, same string as in 2020

# Hits: 40 (Redone in 2024 = 44)

# Limit to Scholarly type, reports
# Hits: 40 (Redone in 2024 = 44)

# NOT mainsubject(freshwater OR bacteria OR aquaculture OR lakes OR rivers OR cyanobacteria OR biogeochemistry)
# Hits: 29 (Redone in 2024 = 38)

ASFA2022 <- synthesisr::read_ref('Data/Bib_files/Input/asfa_search_Jan2022.ris')

# Combine 2022 titles
combo_WoS_ASFA2022 <- bind_rows(WoS2022, ASFA2022)
# Hits: 55

# 2) remove duplicates

# remove duplicates from this combined list of both 2022 searches
combo_match2022 <- revtools::find_duplicates(combo_WoS_ASFA2022, match_variable = "title", group_variables = NULL, match_function = "fuzzdist", method = "fuzz_m_ratio", threshold = 0.1, to_lower = TRUE, remove_punctuation = FALSE)

combo_unique2022 <- extract_unique_references(combo_WoS_ASFA2022, combo_match2022)
# 51 hits (Actually 50 hits????)
# this should be just the new papers published between August 2020 and Jan 2022


# 3) Screen 2020- 2022 titles and abstracts using 'revtools' GUI:
# combo_title_screened2022 <- screen_titles(combo_unique2022) 

#make the screened file into a bib
#write_bibliography(combo_title_screened2022, 'Data/Bib_files/Output/combo_titles_screened_2022.bib', format = 'bib')

#import that bib
all_title_screened2022 <- read_bibliography('Data/Bib_files/Output/combo_titles_screened_2022.bib')

#remove the excluded titles
combo_subset2022 <- all_title_screened2022 %>%
  filter(screened_titles != "excluded")
# Screening the titles --> 18 hits

# Screen the abstracts for inclusion, same criteria

#abstracts_screened2022 <- screen_abstracts(combo_subset2022)

#write_bibliography(abstracts_screened2022, 'Data/Bib_files/Output/abstracts_screened_2022.bib', format = 'bib')
abstracts_screened2022 <- read_bibliography("Data/Bib_files/Output/abstracts_screened_2022.bib")

abstracts_screened_subset2022 <- abstracts_screened2022 %>%
  filter(screened_abstracts != "excluded")
# Screening abstracts --> 2 hits

#write.csv(abstracts_screened_subset2022, 'Data/Bib_files/Output/final_abstracts2022.csv')

# Add these two papers to the main .xls file!!!!!


# Redo the 2020-2022 search -----

# 1) load bib files
wos2022_redo <- synthesisr::read_ref('Data/Bib_files/Input/2024/wos_search_Jan2022_redone.ris')

asfa2022_redo <- synthesisr::read_ref('Data/Bib_files/Input/2024/asfa_search_Jan2022_redone.ris')

# Combine 2022 titles
combo_wos_asfa2022_redo <- bind_rows(wos2022_redo, asfa2022_redo)

# 2) remove duplicates from this combined list of both 2022 searches
combo_match2022_redo <- revtools::find_duplicates(combo_wos_asfa2022_redo, match_variable = "title", group_variables = NULL, match_function = "fuzzdist", method = "fuzz_m_ratio", threshold = 0.1, to_lower = TRUE, remove_punctuation = FALSE)

combo_unique2022_redo <- revtools::extract_unique_references(combo_wos_asfa2022_redo, combo_match2022_redo)
# 124 hits

# OK and NOW let's see how many "new" hits redoing the search got

# take out italics characters bc those seem to mess things up
combo_unique2022_redo_fixed <- combo_unique2022_redo %>%
  mutate(title = gsub('<i>', '', title),
         title = gsub('</i>', '', title))

combo_unique2022_fixed <- combo_unique2022 %>%
  mutate(title = gsub('<i>', '', title),
         title = gsub('</i>', '', title))


redo_2022 <- combo_unique2022_redo_fixed %>%
  filter(!(title %in% combo_unique2022_fixed$title)) # 78 hits to screen??

# 3) Screen 2020- 2022 titles and abstracts using 'revtools' GUI:
# combo_title_screened2022_redo <- screen_titles(redo_2022) 

#make the screened file into a bib
# write_bibliography(combo_title_screened2022_redo, 'Data/Bib_files/Output/combo_title_screened2022_redo.bib', format = 'bib')

#import that bib
all_title_screened2022_redo <- read_refs('Data/Bib_files/Output/combo_title_screened2022_redo.bib')

#remove the excluded titles
combo_subset2022_redo <- all_title_screened2022_redo %>%
  filter(screened_titles != "excluded") %>%
  filter(title != "Impacts of cage fish culture on environment of the Caspian Sea, Caspian zone") %>%
  mutate(label = 1:n())

# Screening the titles --> 5 hits

# Screen the abstracts for inclusion, same criteria

# abstracts_screened2022_redo <- screen_abstracts(combo_subset2022_redo)

#write_bibliography(abstracts_screened2022_redo, 'Data/Bib_files/Output/abstracts_screened2022_redo.bib', format = 'bib')

abstracts_screened2022_redo <- read_bibliography("Data/Bib_files/Output/abstracts_screened2022_redo.bib")

abstracts_screened_subset2022_redo <- abstracts_screened2022_redo %>%
  filter(screened_abstracts != "excluded")
# Screening abstracts -->  2 hits that actually look like real hits..... 


#write_csv(abstracts_screened_subset2022_redo, 'Data/Bib_files/Output/final_abstracts2022_redo.csv')


# 2022/01/04 - 2024/08/27 ------

# 1) Load bib files

# Web of Science search on August 27, 2024

# 303 hits
# Just articles nd proceedings = 280 hits
# Just the categories = 219 hits

wos2024 <- synthesisr::read_ref("Data/Bib_files/Input/2024/wos_search_Aug2024.ris")

# ASFA seach on August 27, 2024

# 42 hits
# Scholarly journals and reports = 42 hits
# exclude categories aquaculture, bacteria, biogeochemistry, cyanobacteria, freshwater, lakes, rivers
# Hits = 23

asfa2024 <- synthesisr::read_ref("Data/Bib_files/Input/2024/asfa_search_Aug2024.ris")

# Combine 2024 titles
combo_wos_asfa2024 <- bind_rows(wos2024, asfa2024)
# Hits: 242

# 2) Remove duplicates 

# remove duplicates from this combined list 
combo_match2024 <- revtools::find_duplicates(combo_wos_asfa2024, match_variable = "title", group_variables = NULL, match_function = "fuzzdist", method = "fuzz_m_ratio", threshold = 0.1, to_lower = TRUE, remove_punctuation = FALSE)

combo_unique2024 <- revtools::extract_unique_references(combo_wos_asfa2024, combo_match2024)
# 227 hits
# this should be just the new papers published between Jan 2022 and August 2024


# 3) Screen 2022- 2024 titles and abstracts using 'revtools' GUI:

# Criteria for inclusion/exclusion:
# [1] Marine environment                             Y/N, yes included
# [2] Animals are the source of nutrients?           Y/N, yes included
# [3] Primary producer measured in response          Y/N, yes included

combo_title_screened2024b <- screen_titles(combo_unique2024) 

#make the screened file into a bib
#synthesisr::write_ref(combo_title_screened2024, 'Data/Bib_files/Output/combo_titles_screened_2024.bib', format = 'bib')

#import that bib
all_title_screened2024 <- synthesisr::read_refs('Data/Bib_files/Output/combo_titles_screened_2024.bib')

#remove the excluded titles
combo_subset2024 <- all_title_screened2024 %>%
  filter(screened_titles != "excluded") %>%
  mutate(screened_abstracts = NA,
         notes = NA,
         label = 1:n()) %>%
  droplevels()
# Screening the titles --> 47 hits

# Screen the abstracts for inclusion, same criteria

abstracts_screened2024 <- revtools::screen_abstracts(combo_subset2024)

#write_bibliography(abstracts_screened2024, 'Data/Bib_files/Output/abstracts_screened_2024.bib', format = 'bib')

# load that file
abstracts_screened2024 <- synthesisr::read_ref("Data/Bib_files/Output/abstracts_screened_2024.bib")

abstracts_screened_subset2024 <- abstracts_screened2024 %>%
  filter(screened_abstracts != "excluded")
# Screening abstracts --> 17 hits! wow!

# write_csv(abstracts_screened_subset2024, 'Data/Bib_files/Output/final_abstracts2024.csv')





# When rerunning the initial searches I get more hits from the 1900-2020 and 2020-2022 searches than I should. I think I should do the 2022-2024 search "normally" and then redo the 1900-2020 and 2020-2022 searches and screen for duplicates and see if there's any usable papers in there and then update the original search numbers to reflect the redone 1900-2020 and 2020-2022 searches

# Rejected papers -----
reject <- read_csv("Data/Bib_files/Output/rejected_papers.csv") %>%
  distinct() # 79 rejected papers

reject_count <- reject%>%
  count(reason) # 5 were duplicates.... not sure where to put those
  


# Graveyard ------
# Exploring why I got diff numbers of papers for redone searches -----

# bind the 2024 wos searches with the original 2020 search
wos_2020_combo_redo <- bind_rows(wos_redo, WoS) 

# which files from this redone search weren't caught in the 2020 search
match_2020_redo <- revtools::find_duplicates(wos_2020_combo_redo, match_variable = "title", group_variables = NULL, match_function = "fuzzdist", method = "fuzz_m_ratio", threshold = 0.1, to_lower = TRUE, remove_punctuation = FALSE)

# create list of just the unique hits between these two searches
combo_unique_2020 <- extract_unique_references(wos_2020_combo_redo, match_2020_redo)

# Which papers did I get when I redid the search vs what did I get in 2020
just_redo_2020 <- combo_unique_2020 %>%
  filter(!(title %in% WoS$title))

# HMMMMM. 
# In 2020, I got 1006 hits, when I redid that search in 2024 I got 1067
# But when I combine those and remove duplicates, I get 1380
# When I just take the refs from the redone search that weren't in the original search, I get 376 "new" refs. 376 + 1006 = 1382, which is almost the 1380 I got when I combined the searches and removed duplicates
# One option is that the titles aren't getting fully screened out
# The other is that papers from the original search are no longer available, and the redone search found 376 new ones, so I have more of a discrepency than the initial 1006 vs 10667

# Now I need to rescreen those 376 papers!!!!! 

# What's happening with ASFA
# bind the asfa 2024 searches with the original 2020 search
asfa_2020_combo_redo <- bind_rows(asfa_redo, ASFA) 

# now see which files from this redone search weren't caught in the 2020 search
asfa_match_2020_redo <- revtools::find_duplicates(asfa_2020_combo_redo, match_variable = "title", group_variables = NULL, match_function = "fuzzdist", method = "fuzz_m_ratio", threshold = 0.1, to_lower = TRUE, remove_punctuation = FALSE)

asfa_combo_unique_2020 <- extract_unique_references(asfa_2020_combo_redo, asfa_match_2020_redo)
# this only results in 718 papers...???

# Which papers did I get when I redid the search vs what did I get in 2020
just_asfa_redo_2020 <- asfa_combo_unique_2020 %>%
  filter(!(title %in% ASFA$title))
# Ok that's just 50 "new ASFA" papers
