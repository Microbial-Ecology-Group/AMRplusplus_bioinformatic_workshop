# Welcome to the "staging" script for statistical analysis of bioinformatic results

# From this script, you can go over what the MEG lab considers to be the main steps for data exploration
# and statistical analysis.
# Remember, the analysis will always have to be based on your study design and performed with the
# goal of testing apriori hypotheses. These scripts are merely meant to provide a guideline.

# The main steps of data analysis we will cover are:
# 1) Loading count matrix results from bioinformatic analyses into R
# 2) Calculating summary statistics
# 3) Creating exploratory figures
# 4) Running some common statistical tests



########################################################################
#                                                                      #
## 1) Loading count matrix results from bioinformatic analyses into R  #
#                                                                      #
########################################################################

# The details of how to load your results into R will depend on how you performed the bioinformatic analysis
# Here, we'll go over how to load resistome and microbiome results that you get from analyzing 
# shotgun metagenomic sequencing reads using the AMR++ pipeline.

# In addition, we'll show you how to load microbiome results from analysis of 16S rRNA sequencing
# reads with qiime2.

# First, we'll need to load the sample "metadata" file. This file must contain the names for all
# of your samples and any other relevant information to help in our analysis.
sample_metadata <- read.table('./data/sample_metadata.csv', header=T, sep=',', row.names = 1, quote = "")
sample_metadata # Run the object name to get more information about the file we just loaded

# Change name of the first column
colnames(sample_metadata)[1] <- "Sample_name"

# Now, let's load the kraken microbiome results.
# The command below runs all of the R code in the script, "scripts/load_kraken_microbiome_data.R"
source("scripts/load_kraken_microbiome_data.R") 
kraken_microbiome.ps # This phyloseq object contains the microbiome results, taxa table, and metadata

# Next, we'll perform similar steps for the resistome data with another script
source("scripts/load_megares_resistome_data.R")
amr.ps # This phyloseq object contains the resistome results, taxa table, and metadata

# Finally, we'll do the same for the qiime2 microbiome results
source("scripts/load_qiime2_microbiome_data.R")
qiime_microbiome.ps # This phyloseq object contains the microbiome results, taxa table, and metadata

#
## To see the details of how to load your data, explore each of the scripts to identify the differences.
#
# We now have these 3 objects to use for further analysis:
# kraken_microbiome.ps
# amr.ps
# qiime_microbiome.ps


########################################################################
#                                                                      #
## 2) Calculating summary statistics                                   #
#                                                                      #
########################################################################








########################################################################
#                                                                      #
## 3) Creating exploratory figures                                     #
#                                                                      #
########################################################################






########################################################################
#                                                                      #
## 4) Running some common statistical tests                            #
#                                                                      #
########################################################################
