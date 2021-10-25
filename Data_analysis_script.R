## Start with this staging file to set up your analysis. 
# Source the utility functions file, which should be in the scripts folder with this file
source('scripts/meg_utility_functions.R')
source('scripts/load_libraries.R')

# Set working directory to the MEG_R_metagenomic_analysis folder and add your data to that folder
#setwd("")

# Set the output directory for graphs:
graph_output_dir = 'graphs'
# Set the output directory for statistics:
stats_output_dir = 'stats'
# In which column of the metadata file are the sample IDs stored?
sample_column_id = 'ID'

####################
## File locations ##
####################
## The files you want to use for input to this (for the MEG group analyses)
## is the AMR_analytic_matrix.csv. So you should have pulled these files from the output of the nextflow pipeline
## and you are now performing this analysis on your local machine. 

## For the AMR analysis, you will also need to download the megares_annotations.csv
## file from the MEGARes website; the annotation file must be from the same version
## of the database as the file you used in the AmrPlusPlus pipeline, i.e. the headers
## must match between the annotation file and the database file.

# Where is the metadata file stored on your machine?
amr_metadata_filepath = 'Proj_3_4_combined_metadata_full.csv'
#amr_count_matrix_filepath = 'megv2_projs_3_4_AMR_analytic_matrix.csv'
amr_count_matrix_filepath = '05AUG21.dedup_megv2_drugs_woSNP_AMR_analytic_matrix.csv' #dedup_megv2_drugs_woSNP_
# Name of the megares annotation file used for this project
megares_annotation_filename = 'data/amr/megares_full_annotations_v2.0.csv'

#################################
## Microbiome - 16S or kraken? ##
#################################

# Where is the metadata file for the microbiome samples stored on your machine?
microbiome_temp_metadata_file <- "Proj_3_4_combined_metadata_full.csv"
# Now, specify file location for 16S
biom_file <- "16S_analysis/exported-biom-table/otu_table_json.biom"
tre_file <- "16S_analysis/exported-tree/tree.nwk"
tax_fasta <- "16S_analysis/exported-rep-seqs/dna-sequences.fasta" #https://data.qiime2.org/2017.6/tutorials/training-feature-classifiers/85_otus.fasta
taxa_file <- "16S_analysis/exported-biom-table-taxa/taxonomy.tsv" #https://data.qiime2.org/2017.6/tutorials/training-feature-classifiers/85_otu_taxonomy.txt


## First, the 16S files. 
# These are the files you'll need to export from qiime2
#qiime tools export project-taxonomy.qza --output-dir exported-biom-table-taxa
#qiime tools export project-rep-seqs.qza --output-dir exported-rep-seqs
#qiime tools export project-aligned-masked-rooted.qza --output-dir exported-tree
#qiime tools export project-dada-table-filtered.qza --output-dir exported-biom-table
#then you need to convert the biom file to "json" using qiime1
#biom convert -i feature-table.biom -o otu_table_json.biom --table-type="OTU table" --to-json




###################
## User Controls ##
###################
## Hopefully, this section should be the only code you need to modify.
## However, you can look into the code in further sections if you need
## to change other, more subtle variables in the exploratory or
## statistical functions.

# The following is a list of analyses based on variables in 
# your metadata.csv file that you want
# to use for EXPLORATORY analysis (NMDS, PCA, alpha rarefaction, barplots)
# NOTE: Exploratory variables cannot be numeric. 

AMR_exploratory_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'Time',
    subsets = list(),
    exploratory_var = 'Time',
    order= c("1st time point","2nd time point", "3rd time point")
  ), 
  list(
    name = 'Time_Group',
    subsets = list(),
    exploratory_var = 'Group',
    order= c("Individual - 1st time point","Individual - 2nd time point","Pen - 1st time point","Pen - 2nd time point","Pen - 3rd time point")
  ),
  # Analysis 7
  # Description: 
  list(
    name = 'Feedlot_cat',
    subsets = list(),
    exploratory_var = 'Feedlot_cat',
    order= ''
  ),
  # Analysis 4
  # Description: 
  list(
    name = '2nd time point_DOF_categories',
    subsets = list('Time != 1st time point'),
    exploratory_var = 'DOF_categories',
    order= c("DOF_1","DOF_2","DOF_3","DOF_4","DOF_5")
  ),
  # Analysis 7
  # Description: 
  list(
    name = '2nd time point_Feedlot_categories',
    subsets = list('Time != 1st time point'),
    exploratory_var = 'Feedlot_categories',
    order= ''
  ),
  list(
    name = 'proj3_2nd time Total_ADD_category',
    subsets = list('Time != 1st time point', 'Project == Project_3'),
    exploratory_var = 'Total_ADD_category',
    order= c("Low total ADD exposure","Medium total ADD exposure","High total ADD exposure")
  ),
  list(
    name = 'proj4_2nd time Total_ADD_category',
    subsets = list('Time != 1st time point', 'Project == Project_4'),
    exploratory_var = 'Total_ADD_category',
    order= c("Low total ADD exposure","Medium total ADD exposure","High total ADD exposure")
  )
)

microbiome_exploratory_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'Time',
    subsets = list(),
    exploratory_var = 'Time',
    order= c("1st time point","2nd time point", "3rd time point")
  ), 
  list(
    name = 'Time_Group',
    subsets = list(),
    exploratory_var = 'Group',
    order= c("Individual - 1st time point","Individual - 2nd time point","Pen - 1st time point","Pen - 2nd time point","Pen - 3rd time point")
  ),
  # Analysis 7
  # Description: 
  list(
    name = 'Feedlot_cat',
    subsets = list(),
    exploratory_var = 'Feedlot_cat',
    order= ''
  ),
  # Analysis 4
  # Description: 
  list(
    name = '2nd time point_DOF_categories',
    subsets = list('Time != 1st time point'),
    exploratory_var = 'DOF_categories',
    order= c("DOF_1","DOF_2","DOF_3","DOF_4","DOF_5")
  ),
  # Analysis 7
  # Description: 
  list(
    name = '2nd time point_Feedlot_categories',
    subsets = list('Time != 1st time point'),
    exploratory_var = 'Feedlot_categories',
    order= ''
  ),
  list(
    name = 'proj3_2nd time Total_ADD_category',
    subsets = list('Time != 1st time point', 'Project == Project_3'),
    exploratory_var = 'Total_ADD_category',
    order= c("Low total ADD exposure","Medium total ADD exposure","High total ADD exposure")
  ),
  list(
    name = 'proj4_2nd time Total_ADD_category',
    subsets = list('Time != 1st time point', 'Project == Project_4'),
    exploratory_var = 'Total_ADD_category',
    order= c("Low total ADD exposure","Medium total ADD exposure","High total ADD exposure")
  )
)

# Each analyses you wish to perform should have its own list in the following
# statistical_analyses list.  A template is provided to get you started.
# Multiple analyses, subsets, and contrasts are valid, but only one random
# effect can be used per analysis.  The contrasts of interest must have their
# parent variable in the model matrix equation.  Contrasts are named by
# parent variable then child variable without a space inbetween, for example:
# PVar1Cvar1 where the model matrix equation is ~ 0 + Pvar1.
AMR_statistical_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'proj3_over_time',
    subsets = list('Project == Project_3'),
    model_matrix = '~ 0 + Time',
    contrasts = list('Time2nd time point-Time1st time point'),
    random_effect = NA
  ),
  # Analysis 2
  # Description: 
  list(
    name = 'proj4_over_time',
    subsets = list('Project == Project_4', 'Time != 3rd time point'),
    model_matrix = '~ 0 + Time',
    contrasts = list('Time2nd time point-Time1st time point'),
    random_effect = NA
  )
)

microbiome_statistical_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'proj3_over_time',
    subsets = list('Project == Project_3'),
    model_matrix = '~ 0 + Time',
    contrasts = list('Time2nd time point-Time1st time point'),
    random_effect = NA
  ),
  # Analysis 2
  # Description: 
  list(
    name = 'proj4_over_time',
    subsets = list('Project == Project_4', 'Time != 3rd time point'),
    model_matrix = '~ 0 + Time',
    contrasts = list('Time2nd time point-Time1st time point'),
    random_effect = NA
  )
)



### Run the analysis
##
#
## Pick the correct script that handles resistome data and/or microbiome data. 
# Choose one of the following scripts
source('scripts/metagenomeSeq_qiime2.R')
source('scripts/metagenomeSeq_megaresv2_onlydrugs.R')


# Metadata
source('scripts/1_Metadata_editing.R')
# Add metadata to melted tables
setkey(amr_melted_raw_analytic,ID) 
setkey(amr_melted_analytic,ID) 
setkey(microbiome_melted_analytic,ID)
setkey(microbiome_melted_raw_analytic,ID)

setkey(metadata,ID)
setkey(microbiome_metadata,ID)

amr_melted_analytic <- amr_melted_analytic[metadata]
amr_melted_raw_analytic <- amr_melted_raw_analytic[metadata]
microbiome_melted_analytic <- microbiome_melted_analytic[microbiome_metadata]
microbiome_melted_raw_analytic <- microbiome_melted_raw_analytic[microbiome_metadata]

full_metadata <- cbind(microbiome_metadata, metadata)
write.csv(full_metadata,"updated_full_metadata.csv")

write.csv(microbiome_metadata,"diversity_microbiome_metadata.csv")
write.csv(metadata,"diversity_metadata.csv")

# After running this script, these are the useful objects that contain all the data aggregated to different levels
# The metagenomeSeq objects are contained in these lists "AMR_analytic_data" and "microbiome_analytic_data"
# Melted counts are contained in these data.table objects "amr_melted_analytic" "microbiome_melted_analytic"

# Run code to make some exploratory figures, zero inflated gaussian model, and output count matrices.
## Choose one of the following scripts:
source('scripts/print_microbiome_figures.R')
source('scripts/print_AMR_figures.R')

## For ZIG model results
#source('scripts/print_microbiome_ZIG_results.R')
source('scripts/print_AMR_ZIG_results.R')


