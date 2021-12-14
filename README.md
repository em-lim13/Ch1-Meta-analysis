# Ch1-Meta-analysis

This repository stores all the code and data used for the manuscript "Consumer-mediated nutrient cycling in marine environments: a global meta-analysis," accepted by Global Ecology and Biogeography (October 2026). The full paper is available here:


## Code
Lit_Search_Script.R: 
This file takes the results of the Web of Science and ASFA literature searches (Data/Bib_files/Input) and combines them, removes duplicates, and then uses revtools to screen the titles and abstracts for inclusion. At each stage of screening I write .bib files to Bib_files/Output so I can jump back into the process without having to rerun everything. The final output is the final_titles.csv, which I used as the base csv for the data extraction phase

Mean_se_calc.R: 
This file uses the data in the Data/Mean_se_calc_data folder to calculate means and standard error for a subset of publications using their published data.

Functions.R: 
This file stores all of the functions I wrote to avoid copy and pasting the same things again and again. Plotting themes are also stored here. 

SST_calc.R:
This script calculates downloads sea surface temperature data for each study location and duration and writes this data to the file "meta_sst_final.csv" 

LME_calc.R:
This script uses the large marine ecosystem polygons stored in the "lme_shp" folder to determine which region each of the studies was conducted in, and writes this data to the file "LME_df.csv"

Meta-analysis_Script.R: 
This is the main code file! This uses the meta_analysis_data.csv created in Lit_Search_Script.R where all of the meta-data is stored. It also sources the functions stored in Functions.R. This is where we manipulate the data, calculate effect sizes, create and analyse the models, create figures and tables, and calculate summary statistics

## Data:
meta_data.docx:
This is the meta data for meta_analysis_data2025.csv where I explain all of the variables

meta_analysis_data.csv: This is the main csv for the meta-analysis

meta_analysis_data.xlsx: This is the raw excel file for the meta-analysis

Bib_files: folder of files generated and used in Lit_Search_Script.R

lme_shp: folder of large marine ecoregion shapefiles

Mean_SE_calc_data: raw data from publications included in the meta-analysis to calculate mean and standard error

Output_data: folder where .csv's written in scripts for use in other scripts live

## Figures
All figures generated for the publication

## Tables 
All tables generated for the publication

