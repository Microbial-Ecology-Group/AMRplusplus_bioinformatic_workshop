# Welcome to the "staging" script for statistical analysis of bioinformatic results

# From this script, you can go over what the MEG lab considers to be the main steps for data exploration
# and statistical analysis.
# Remember, the analysis will always have to be based on your study design and performed with the
# goal of testing apriori hypotheses. These scripts are merely meant to provide a guideline.

# The main steps of data analysis we will cover are:
# 1) Loading count matrix results from bioinformatic analyses into R
# 2) Calculating summary statistics
# 3) Normalizing counts and creating exploratory figures
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
# Here's an example of the file we are loading in:
kraken_microbiome <- read.table('./data/kraken_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")
# Notice samples are on the columns and taxa are the rows
colnames(kraken_microbiome)
rownames(kraken_microbiome)

# The command below runs all of the R code in the script, "scripts/load_kraken_microbiome_data.R"
source("scripts/load_kraken_microbiome_data.R") 
kraken_microbiome.ps # This phyloseq object contains the microbiome results, taxa table, and metadata

# Next, we'll perform similar steps for the resistome data with another script
source("scripts/load_megares_resistome_data.R")
amr.ps # This phyloseq object contains the resistome results, taxa table, and metadata

# Finally, we'll do the same for the qiime2 microbiome results
source("scripts/load_qiime2_microbiome_data.R")
qiime_microbiome.ps # This phyloseq object contains the microbiome results, taxa table, and metadata


# We now have these 3 phyloseq objects to use for further analysis:
# kraken_microbiome.ps
# amr.ps
# qiime_microbiome.ps

# You can access the various components of the phyloseq object using the following functions
sample_data(kraken_microbiome.ps)
otu_table(kraken_microbiome.ps)
tax_table(kraken_microbiome.ps)

#
## To see the details of how to load your data, explore each of the scripts to identify the differences.
#

########################################################################
#                                                                      #
## 2) Calculating summary statistics                                   #
#                                                                      #
########################################################################

# For steps 2-4, we'll use the kraken microbiome results to go over various examples
# If you want a quick way to replicate this analysis for the other data types, highlight everything below
# and change the "kraken_microbiome.ps" object to either "amr.ps" for the resistome data or 
# "qiime_microbiome.ps" for the microbiome data.

# This is using columns in your sample metadata file, so be sure it corresponds with your data
# Summarize for entire kraken microbiome dataset

# Notice we can access the metadata file columns by combining the "sample_data()" function and the $
# to call out a specific column, in this case "X16S_Raw_paired_reads"
sample_data(kraken_microbiome.ps)$X16S_Raw_paired_reads

# We can make simple calculations using built-in R functions like "sum(), mean(), etc"
# This is the total sum of raw paired reads across our entire dataset
sum(sample_data(kraken_microbiome.ps)$X16S_Raw_paired_reads)

# Next, we might want to make multiple calculations on our metadata variables
# Here, we can use the "dplyr" package with the "%>%" function and the "summarize()" function as shown below:
sample_data(kraken_microbiome.ps) %>%
  summarize(total_16S_counts = sum(X16S_Raw_paired_reads), mean_16S_counts = mean(X16S_Raw_paired_reads),
            min_16S_counts = min(X16S_Raw_paired_reads),max_16S_counts = max(X16S_Raw_paired_reads),
            median_16S_counts = median(X16S_Raw_paired_reads))

# We can group our summary statistics by specific metadata variables by adding the "group_by()" function
# Ignore the warning about the R class, this shouldn't affect your downstream analysis
sample_data(kraken_microbiome.ps) %>%
  group_by(Group) %>% 
  summarize(total_16S_counts = sum(X16S_Raw_paired_reads), mean_16S_counts = mean(X16S_Raw_paired_reads),
            min_16S_counts = min(X16S_Raw_paired_reads),max_16S_counts = max(X16S_Raw_paired_reads),
            median_16S_counts = median(X16S_Raw_paired_reads))

# We have many other metadata variables to use for further analysis such as the # of reads after QC filtering
# total non-host reads, and the number of mapped reads, to name a few. If you write up your data for 
# submission to a peer-reviewed journal, you'll often describe these results in the text and provide a
# table of the raw values.

#
## Agglomerate ASV counts to different taxonomic levels
#

# Currently, our microbiome counts are not aggregated to any taxonomic levels. This means that one row
# can be labeled down to the species level while another can only go to the phylum level.


# Using tax_glom(), we can easily aggregate counts to different levels (taxonomic or AMR annotation levels)
kraken_phylum.ps <- tax_glom(kraken_microbiome.ps, "phylum")
kraken_phylum.ps

#
## To get more ideas and explore further ways to summarize your data, 
# explore the script, "scripts/Step2_R_Data_summary_statistics.R"
#

########################################################################
#                                                                      #
## 3) Normalizing counts and creating exploratory figures              #
#                                                                      #
########################################################################

# As you might have noticed with the summary statistics, there was variation in both the 
# number of sequenced and mapped reads between samples. This could be an effect of sequencing bias and
# you must account for this using count normalization. The MEG lab currently uses cumulative sum scaling (CSS)

# In this section, we'll go over how to normalize the microbiome results and create some exploratory
# figures to further explore your data.

# Just visually, we can observe differences in the number of total mapped reads between samples
plot_raw_kraken_phylum <- plot_bar(kraken_phylum.ps, fill = "phylum") + 
  facet_wrap(~ Group, scales = "free_x") +
  labs(title= "Raw kraken microbiome counts") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45),
        panel.background = element_blank(),
        legend.title = element_text( size=8), 
        legend.text=element_text(size=8))
plot_raw_kraken_phylum

# Notice we can change the "facet_wrap" to group our data by phyla, instead of group
plot_raw_kraken_phylum_counts <- plot_bar(kraken_phylum.ps, fill = "phylum") + 
  facet_wrap(~ phylum, scales = "free_x") +
  labs(title= "Raw kraken microbiome counts by phyla") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45),
        panel.background = element_blank(),
        legend.title = element_text( size=8), 
        legend.text=element_text(size=8))
plot_raw_kraken_phylum_counts

# Prior to using the full dataset for count normalization, we might choose a threshold of 
# counts that each taxa much meet. For a simple example, we'll use "50" counts as the minimum threshold.
filtered_microbiome.ps = filter_taxa(kraken_microbiome.ps, function(x) sum(x) > 50, TRUE)
# First, we convert the phyloseq object to metagenomeSeq
library(metagenomeSeq)
kraken_microbiome.metaseq <- phyloseq_to_metagenomeSeq(filtered_microbiome.ps)

# Check out this object:
kraken_microbiome.metaseq

# This is where we perform the normalization
kraken_microbiome.metaseq <- cumNorm(kraken_microbiome.metaseq)

# Like phyloseq, metagenomeSeq has it's own functions for accessing their data.
# Here, we need to use MRcounts() and re-make the phyloseq object with the normalized counts
CSS_kraken_microbiome_counts <- MRcounts(kraken_microbiome.metaseq, norm = TRUE)

# Use the new counts and merge with components from our original phyloseq object.
CSS_normalized_kraken.ps <- merge_phyloseq(otu_table(CSS_kraken_microbiome_counts, taxa_are_rows = TRUE),sample_data(filtered_microbiome.ps),tax_table(filtered_microbiome.ps), phy_tree(filtered_microbiome.ps))

# Aggregate counts to phylum
CSS_normalized_phylum_qiime.ps <- tax_glom(CSS_normalized_qiime.ps, "phylum")

# Notice the y-axis values are counts and not proportions.
plot_css_qiime_phylum <- plot_bar(CSS_normalized_phylum_qiime.ps, fill = "phylum") + 
  facet_wrap(~ Group, scales = "free_x") +
  labs(title= "CSS normalized qiime microbiome counts") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45),
        panel.background = element_blank(),
        legend.title = element_text( size=8), 
        legend.text=element_text(size=8))
plot_css_qiime_phylum




########################################################################
#                                                                      #
## 4) Running some common statistical tests                            #
#                                                                      #
########################################################################

# The statistical analysis portion will be the most varied and 
