if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")

# Following : https://huttenhower.sph.harvard.edu/maaslin
#Load MaAsLin2 package into the R environment
library(Maaslin2)
?Maaslin2

rm(list=ls())
getwd()
setwd("~/marine_drugs/marine_drugs/reports/figures/")
dir.create("pathway-profiling") # Create a new directory
setwd("pathway-profiling") # Change the current working directory 
getwd() #check if directory has been successfully changed


# MaAsLin2 requires two input files
# one for taxonomic or functional feature abundances
# and one for sample metadata.

### Example multi-variate association analysis using Maaslin2
# 0.1 example data filepaths:
input_data = system.file(
  "extdata", "HMP2_taxonomy.tsv", package="Maaslin2") # The abundance table file
input_data
input_metadata = system.file(
  "extdata", "HMP2_metadata.tsv", package="Maaslin2") # The metadata table file
input_metadata
# 0.2 read in example data filepaths as dataframes
df_input_data = read.table(file = input_data, header = TRUE, sep = "\t",
                          row.names = 1,
                          stringsAsFactors = FALSE)
df_input_data[1:5, 1:5]
df_input_metadata = read.table(file = input_metadata, header = TRUE, sep = "\t",
                               row.names = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, ]
# Run a multivariable regression model to test for 
# the association between microbial species abundance versus IBD diagnosis
fit_data = Maaslin2(
  input_data = input_data, 
  input_metadata = input_metadata, 
  output = "demo_output", 
  fixed_effects = c("diagnosis", "dysbiosis"))

## END OF EXAMPLE ANALYSIS


####
# BEGIN Sponge activity Analysis
###

getwd()

# 1. Data file of pathway abundances
# The abundance table file
# NOTE: We removed the comment character'#' from the header column s.t. read.table(...) appropriately parses the cols.
input_data = "/Users/rees/marine_drugs/marine_drugs/data/processed/pathway-profiling/master_community_genefamilies_profiles.munged.tsv"
input_data

# 2. Metadata file of sponge phenotypes
input_metadata = "/Users/rees/marine_drugs/marine_drugs/data/raw/sponge_metadata.rnaseqsponges.munged.tsv"
input_metadata

# Reading inputs as dataframes
df_data = read.table(file = input_data,
                     header = T,
                     quote = "",
                     row.names = 1,
                     sep = "\t",
                     stringsAsFactors = FALSE,
                     fill=TRUE)
# Preview dataframe
df_data[1:5, 1:5]
# Now read in sponge metadata (bioactivity and taxonomy)
df_metadata = read.table(file = input_metadata,
                         header = TRUE,
                         row.names = 1,
                         sep = "\t",
                         quote = "",
                         stringsAsFactors = FALSE)
# Preview dataframe
df_metadata[1:5, ]

# ?subset
# df_metadata <- subset(df_metadata, Metatranscriptome == "Yes")
# df_metadata

# Initially received error about some library being corrupt:
# Solution: https://github.com/lme4/lme4/issues/407#issuecomment-273537772

fit_data = Maaslin2(
  input_data = df_data, 
  input_metadata = df_metadata, 
  output = "genefamilies_bioactivity_associations", 
  fixed_effects = c("Activity", "Taxonomic.Classification"),
  cores=2)

## Associations between pathway abundances

input_data = "/Users/rees/marine_drugs/marine_drugs/data/processed/pathway-profiling/master_community_pathabundance_profiles.munged.tsv"
input_data

# 2. Metadata file of sponge phenotypes
input_metadata = "/Users/rees/marine_drugs/marine_drugs/data/raw/sponge_metadata.rnaseqsponges.munged.tsv"
input_metadata

# Reading inputs as dataframes
df_data = read.table(file = input_data,
                     header = T,
                     quote = "",
                     row.names = 1,
                     sep = "\t",
                     stringsAsFactors = FALSE,
                     fill=TRUE)
# Preview dataframe
df_data[1:5, 1:5]
# Now read in sponge metadata (bioactivity and taxonomy)
df_metadata = read.table(file = input_metadata,
                         header = TRUE,
                         row.names = 1,
                         sep = "\t",
                         quote = "",
                         stringsAsFactors = FALSE)
# Preview dataframe
df_metadata[1:5, ]

# ?subset
# df_metadata <- subset(df_metadata, Metatranscriptome == "Yes")
# df_metadata

# Initially received error about some library being corrupt:
# Solution: https://github.com/lme4/lme4/issues/407#issuecomment-273537772

fit_data = Maaslin2(
  input_data = df_data, 
  input_metadata = df_metadata, 
  output = "pathway_abundances_bioactivity_associations", 
  fixed_effects = c("Activity", "Taxonomic.Classification"),
  cores=2)