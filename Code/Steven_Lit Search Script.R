##

# Dreissenid HS Model Lit Review Script
# Date: July 7 2019, rev. Jan 06 2020
# Author: Steven Brownlee
# Email: sbrownle@sfu.ca

# Description: Script using the package 
#'devtools' to manage and process exported
# .bib files retrieved from literature search.

##

setwd("C:\\SFU\\AA - PRISMA Lit Review\\06-01-2020 Search")

# 1.) Install package as needed, load library:

#install.packages('revtools')
#install.packages('tidyverse')
#install.packages('devtools')
#install.packages('wosr')

library(revtools)
library(tidyverse)
library(wosr)

collections <- c('SCI', 'SSCI', 'ISTP', 'ISSHP')

###
sid <- auth()

query <- "(mussel* AND (dreissen* OR zebra OR quagga) AND (invas* OR non-native OR introduc*))
# AND (habitat OR environment* OR niche) AND (suitability OR requirement* OR limit* OR distribut*)"

pull_wos()

###

# 2.) Load exported .bib/.ris file from database.

# Search terms: (mussel* AND (dreissen* OR zebra OR quagga) AND (invas* OR non-native OR introduc*))
# AND (habitat OR environment* OR niche) AND (suitability OR requirement* OR limit*)

# Web of Science search, January 06 2020:

WoS_01 <- read_bibliography('wos_search_2020.bib')

# write.csv(WoS_01, 'WoS_01.csv')

# ASFA search, January 07 2020: 

ASFA_a <- read_bibliography('asfa_search_2020_a.ris')
ASFA_b <- read_bibliography('asfa_search_2020_b.ris')
ASFA_c <- read_bibliography('asfa_search_2020_c.ris')
ASFA_d <- read_bibliography('asfa_search_2020_d.ris')

ASFA_comb <- rbind(ASFA_a, ASFA_b, ASFA_c, ASFA_d)

# write.csv(ASFA_01_comb, 'ASFA_01.csv')

# 3.) Identify and remove duplicates based on title:

WoS_match <- find_duplicates(WoS_01, match_variable = "title", group_variables = NULL, match_function = "fuzzdist",
                             method = "fuzz_m_ratio", threshold = 0.1, to_lower = TRUE, remove_punctuation = FALSE)
ASFA_match <- find_duplicates(ASFA_comb, match_variable = "title", group_variables = NULL, match_function = "fuzzdist",
                              method = "fuzz_m_ratio", threshold = 0.1, to_lower = TRUE, remove_punctuation = FALSE)

WoS_unique1 <- extract_unique_references(WoS_01, WoS_match)
ASFA_unique1 <- extract_unique_references(ASFA_comb, ASFA_match)

# save(WoS_unique_2020, file = "WoS_unique.rData")

# save(ASFA_unique_2020, file = "ASFA_unique.rData")

load("WoS_unique_2020.rData")
load("ASFA_unique_2020.rData")

# 4.) Screen titles and abstracts using 'revtools' GUI:

# Criteria for inclusion/exclusion:

# [1.] Zebra/quagga mussels topic of exploration?     Y/N, yes included
# [2.] Native or non-native range?                    Non-native included
# [3.] Observational or experimental?                 Observational included
# [4.] Does it discuss correlates of distribution?    Y/N, yes included

#

WoS_title_screened <- screen_titles(WoS_unique) 

# write_bibliography(WoS_title_screened, 'WoS_title_screened_2020.bib', format = 'bib')

WoS_title_screened <- read_bibliography('WoS_title_screened_2020.bib')

WoS_subset <- WoS_title_screened[WoS_title_screened$selected != 'excluded',]

WoS_abstract_screened <- screen_abstracts(WoS_subset)

# write_bibliography(WoS_abstract_screened, 'WoS_abstract_screened_2020.bib', format = 'bib')
# WoS_abstract_screened <- read_bibliography("WoS_abstract_screened_2020.bib")

WoS_abstract_screened_2020 <- read_bibliography('WoS_abstract_screened_2020.bib')

WoS_abstract_screened_subset <- WoS_abstract_screened_2020[WoS_abstract_screened_2020$selected != 'excluded',]

WoS_titles <- WoS_abstract_screened_subset %>% select(title)

#

ASFA_title_screened <- screen_titles(ASFA_unique)

# write_bibliography(ASFA_title_screened, 'ASFA_title_screened.bib', format = 'bib')

ASFA_title_screened <- read_bibliography('ASFA_title_screened.bib')

ASFA_subset<- ASFA_title_screened[ASFA_title_screened$screened_titles != 'excluded',]

ASFA_abstract_screened <- screen_abstracts(ASFA_subset)

# write_bibliography(ASFA_abstract_screened, 'ASFA_abstract_screened.bib', format = 'bib')

ASFA_abstract_screened <- read_bibliography('ASFA_abstract_screened_2020.bib')

# # # ASFA_abstract_screened['8', 'selected'] <- 'excluded'

ASFA_abstract_screened_subset <- ASFA_abstract_screened[ASFA_abstract_screened$screened_abstracts != 'excluded',]

ASFA_titles <- ASFA_abstract_screened_subset %>% select(title)

## ASFA_test <- ASFA_abstract_screened_subset[ASFA_abstract_screened_subset$notes == 'Include, ambiguous.',]
## WoS_test <- WoS_abstract_screened_subset[WoS_abstract_screened_subset$notes == 'Include, ambiguous.',]

# 5.) Concatenate final list of papers from ASFA and WoS, filter for duplicates and retrieve 
# records for final review list.

titles_comb <- rbind(WoS_titles, ASFA_titles)

titles_unique <- unique(titles_comb)

## write.csv(titles_unique, 'final_titles.csv')

## papers_final <- merge(x = titles_unique, y = ASFA_abstract_screened_subset, by = 'title', all.y = TRUE)
