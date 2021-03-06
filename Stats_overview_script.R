# Welcome to the "staging" script for statistical analysis of bioinformatic results
# Run this script to install all required R packages, or comment it out if not needed
#source('scripts/install_course_packages.R')

# Once the packages are installed, run the command below to load all the packages
source('scripts/load_R_packages.R')

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
# of your samples in the first column and any other relevant information to help in our analysis.
sample_metadata <- read.table('./data/megarich_sample_metadata.csv', header=T, sep=',', row.names = 1, quote = "")
sample_metadata # Run the object name to get more information about the file we just loaded

# Now, let's load the AMR resistome results.
# Here's an example of the file we are loading in:
amr <- read.table('./data/shotgun_AMR_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")
# Notice samples are on the columns and taxa are the rows
colnames(amr)
rownames(amr)
View(amr)

# The command below runs all of the R code in the script, "scripts/Step1_load_megares_resistome_data.R"
source("scripts/Step1_load_megares_resistome_data.R")
amr.ps # This phyloseq object contains the resistome results, taxa table, and metadata

# Next, we'll perform similar steps for the microbiome data with another script
source("scripts/Step1_load_kraken_microbiome_data.R") 
kraken_microbiome.ps # This phyloseq object contains the microbiome results, taxa table, and metadata

# For the target enriched data, the resistome results will be loaded the same way as above
source("scripts/Step1_load_TE_megares_resistome_data.R")
TE_amr.ps # This phyloseq object contains the resistome results, taxa table, and metadata

# Finally, we'll do the same for the qiime2 microbiome results
source("scripts/Step1_load_qiime2_microbiome_data.R")
qiime_microbiome.ps # This phyloseq object contains the microbiome results, taxa table, and metadata


# We now have these 4 phyloseq objects to use for further analysis:
# kraken_microbiome.ps
# amr.ps
# TE_amr.ps
# qiime_microbiome.ps

# You can access the various components of the phyloseq object using the following functions:
sample_data(qiime_microbiome.ps)
otu_table(qiime_microbiome.ps)
tax_table(qiime_microbiome.ps)

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
# and change the "qiime_microbiome.ps" object to either "amr.ps" for the resistome data or 
# "qiime_microbiome.ps" for the microbiome data.

# This is using columns in your sample metadata file, so be sure it corresponds with your data
# Summarize for entire kraken microbiome dataset

# Notice we can access the metadata file columns by combining the "sample_data()" function and the $
# to call out a specific column, in this case "Raw_paired_reads"
sample_data(qiime_microbiome.ps)$Raw_paired_reads

# We can make simple calculations using built-in R functions like "sum(), mean(), etc"
# This is the total sum of raw paired reads across our dataset
sum(sample_data(qiime_microbiome.ps)$Raw_paired_reads)
mean(sample_data(qiime_microbiome.ps)$Raw_paired_reads)
min(sample_data(qiime_microbiome.ps)$Raw_paired_reads)
max(sample_data(qiime_microbiome.ps)$Raw_paired_reads)

# Remember that we used two different metadata files, one for the 16S dataset analyzed with qiime2 and another
# for the 64 metagenomic sequencing samples from the MEGARich dataset (kraken microbiome, resistome, 
# target enriched resistome)
sample_data(kraken_microbiome.ps)

# We have many other metadata variables to use for further analysis such as the # of reads after QC filtering
# total non-host reads, and the number of mapped reads, to name a few. If you write up your data for 
# submission to a peer-reviewed journal, you'll often describe these results in the text and provide a
# table of the raw values.

#
## Alpha diversity measures
#
# We can calculate a bunch of diversity measures with a single command, estimate_richness()
qiime_microbiome_16S_diversity_values <- estimate_richness(qiime_microbiome.ps)
min(qiime_microbiome_16S_diversity_values$Observed)

#
## Agglomerate ASV counts to different taxonomic levels
#

# Currently, our microbiome counts are not aggregated to any taxonomic levels. This means that one row
# can be labeled down to the species level while another can only go to the phylum level.

# Using tax_glom(), we can easily aggregate counts to different levels (taxonomic or AMR annotation levels)
qiime_phylum.ps <- tax_glom(qiime_microbiome.ps, "phylum")
qiime_phylum.ps

# The "Fisher" diversity index will cause an error due the phylum aggregated data not having any "singletons"
phylum_diversity <- estimate_richness(qiime_phylum.ps)

# You can try again, but this time only ask for some of the diversity indices
phylum_diversity <- estimate_richness(qiime_phylum.ps, measures = c("Observed", "Shannon", "InvSimpson"))

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
# figures to further explore your data. First, we'll start with some figures to visualize summary statistics.

# We'll pick up were we left off with the diversity indices
qiime_microbiome_16S_diversity_values
# We need to add a column with the sample names to help us merge this data to the metadata file and
# add it back to the phyloseq object
qiime_microbiome_16S_diversity_values$Name <- row.names(qiime_microbiome_16S_diversity_values)

# Now, we use left_join() to add the sample_metadata to the object with diversity values
expanded_metadata <- left_join(qiime_sample_metadata,qiime_microbiome_16S_diversity_values, by = "Name")

# Need to provide row.names to the newly formed metadata object
row.names(expanded_metadata) <- expanded_metadata$Name

# Merge updated metadata values back into phyloseq object
qiime_microbiome.ps <- merge_phyloseq(qiime_microbiome, phy_tree(qiime_microbiome_phylo_tree), tax_table(as.matrix(taxa.df)), sample_data(expanded_metadata))

# Notice the metadata file has been updated in the phyloseq object
sample_data(qiime_microbiome.ps)

# Using base graphics
# (Y ~ X)
# Create boxplot of raw paired reads grouped by Sample_type
boxplot(sample_data(qiime_microbiome.ps)$Raw_paired_reads ~ sample_data(qiime_microbiome.ps)$Sample_type)

# Now, plot the phred quality scores
boxplot(sample_data(qiime_microbiome.ps)$Mean_phred_score ~ sample_data(qiime_microbiome.ps)$Sample_type)

# Observed diversity (Richness)
boxplot(sample_data(qiime_microbiome.ps)$Observed ~ sample_data(qiime_microbiome.ps)$Sample_type)

#
## Use ggplot to create boxplot of raw paired reads by Sample_type
#
# Notice that we have to specify that the Y value is numeric with "as.numeric()"
ggplot(sample_data(qiime_microbiome.ps), aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_point()

# Like ggplot object, you can keep adding "layers" to modify the script
ggplot(sample_data(qiime_microbiome.ps), aes(x = Sample_type , y = as.numeric(Raw_paired_reads), color = Sample_type)) + 
  geom_boxplot() +
  geom_point() +
  labs(title = "16S Metagenomic sequencing reads", x = "Sample type", y = "raw paired reads") + 
  theme(axis.text.x = element_text( size = 18),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank())

#
## Diversity indices boxplots
#
ggplot(sample_data(qiime_microbiome.ps), aes(x = Sample_type, y = Observed, color = Sample_type)) +
  geom_boxplot() +
  geom_point() +
  labs(title = "Unique features by treatment group", x = "Sample type", y = "Observed features") + 
  theme_classic()


#
## Count normalization
#
# Just visually, we can observe differences in the number of total mapped reads between samples
#
plot_raw_qiime_phylum <- plot_bar(qiime_phylum.ps, fill = "phylum") + 
  facet_wrap(~ Sample_type, scales = "free_x") +
  labs(title= "Raw 16S microbiome counts") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45),
        panel.background = element_blank(),
        legend.title = element_text( size=8), 
        legend.text=element_text(size=8))
plot_raw_qiime_phylum


# Example of rarefaction
# Set seed for reproducibility
set.seed(42)
# Rarefy to sample with lowest number of reads
rarefied_qiime.ps <- rarefy_even_depth(qiime_microbiome.ps, sample.size = min(sample_sums(qiime_microbiome.ps)))
# Aggregate counts to phylum
rarefied_phylum_qiime.ps <- tax_glom(rarefied_qiime.ps, "phylum")
# Bar plot to visualize results
plot_rarefied_qiime_phylum <- plot_bar(rarefied_phylum_qiime.ps, fill = "phylum") + 
  facet_wrap(~ Sample_type, scales = "free_x") +
  labs(title= "Rarefied qiime microbiome counts") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45),
        panel.background = element_blank(),
        legend.title = element_text( size=8), 
        legend.text=element_text(size=8))
plot_rarefied_qiime_phylum


#
## Normalization using Cumulative sum scaling
# First, we convert the phyloseq object to metagenomeSeq
# load the metagenomeSeq library
library(metagenomeSeq)

# We can use the following function to convert object
qiime_microbiome.metaseq <- phyloseq_to_metagenomeSeq(qiime_microbiome.ps)
# Check out this metagenomeSeq object:
qiime_microbiome.metaseq

# This is where we perform the normalization and add the values to the metagenomeSeq object
qiime_microbiome.metaseq <- cumNorm(qiime_microbiome.metaseq)

# Like phyloseq, metagenomeSeq has it's own functions for accessing their data.
# Here, we need to use MRcounts() to extract the normalized counts.
# We'll then re-make the phyloseq object with the normalized counts
CSS_qiime_microbiome_counts <- MRcounts(qiime_microbiome.metaseq, norm = TRUE)
# Use the new counts and merge with components from our original phyloseq object.
CSS_normalized_qiime.ps <- merge_phyloseq(otu_table(CSS_qiime_microbiome_counts, taxa_are_rows = TRUE),sample_data(qiime_microbiome.ps),tax_table(qiime_microbiome.ps))
# Aggregate counts to phylum
CSS_normalized_phylum_qiime.ps <- tax_glom(CSS_normalized_qiime.ps, "phylum")

# Notice the y-axis values are counts and not proportions.
plot_css_qiime_phylum <- plot_bar(CSS_normalized_phylum_qiime.ps, fill = "phylum") + 
  facet_wrap(~ Sample_type, scales = "free_x") +
  labs(title= "CSS normalized qiime microbiome counts") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45),
        panel.background = element_blank(),
        legend.title = element_text( size=8), 
        legend.text=element_text(size=8))
plot_css_qiime_phylum

#
## Here, we can save both figures to .pdf files and compare the results
# You first need to run the object containg the plot, and then the "ggsave()" function right after that.
plot_raw_qiime_phylum
ggsave("raw_qiime_phylum_counts_by_sample.pdf",width = 30, height = 20, units = "cm")
# Repeat the process for another figure
plot_rarefied_qiime_phylum
ggsave("CSSnormalized_qiime_phylum_counts_by_sample.pdf",width = 30, height = 20, units = "cm")
# Repeat the process for another figure
plot_css_qiime_phylum
ggsave("CSSnormalized_qiime_phylum_counts_by_sample.pdf",width = 30, height = 20, units = "cm")

#
## Relative abundance plots
#
# Instead of showing total counts per sample, we often calculate the relative abundance 
# Convert OTU abundances to relative abundances
rel_phylum_qiime.ps <- transform_sample_counts(CSS_normalized_phylum_qiime.ps, function(x) x / sum(x) )

# We can plot these results, notice the y-axis in these relative abundance plots
plot_bar(rel_phylum_qiime.ps, fill= "phylum") +
  facet_wrap(~ Sample_type, scales = "free_x")


#
## Ordination plots
#
# To compare the microbiome composition between samples, we can calculate the "Bray-Curtis" distance
# between samples using the "ordinate()" function.
ordination_phylum_bray <- ordinate(CSS_normalized_phylum_qiime.ps, method = "NMDS")
# We can then use "plot_ordination()" to plot the distance matrix
# We specify that we want to compare "samples" and color the points by the "Sample_type" metadata variable
plot_ordination(CSS_normalized_phylum_qiime.ps, ordination_phylum_bray, type = "samples",color = "Sample_type")


#
## To try out some other methods of count normalization and other types of figures, 
# explore the script, "scripts/Step3_Count_normalization.R", "scripts/Step3_Introduction_to_plotting", 
# and "scripts/Step3_Ordination.R"
# 

########################################################################
#                                                                      #
## 4) Running some common statistical tests                            #
#                                                                      #
########################################################################

# Running statistical analyses on your microbiome and resistome data is by far the most challenging 
# component and has to be tailored to fit your specific dataset and research objectives.

# A common statistical test we can perform is a non-parametric Wilcoxon test to compare the total number
# of sequenced reads between samples groups. This test could similarly be used to test for differences 
# between the total number of mapped reads to the microbiome and resistome, or comparing 

# Test for differences in sequencing depth by Sample_type
# The format is wilcox.test( Y numeric values ~ X grouping factor)
# Here, we'll compare the number of raw reads between Sample_type.
# Notice, we'll get an error that the "grouping factor must have exactly 2 levels"
wilcox.test(sample_data(CSS_normalized_phylum_qiime.ps)$Raw_paired_reads ~ sample_data(CSS_normalized_phylum_qiime.ps)$Sample_type)

# Wilcoxon tests are useful for comparisons between two groups, but you might have to test between 
# more than two groups, such as sequencing lanes, so we can use generalized linear models, glm():
glm(sample_data(CSS_normalized_phylum_qiime.ps)$Raw_paired_reads ~ sample_data(CSS_normalized_phylum_qiime.ps)$Sample_type)
# We need to use "summary()" to view the results
summary(glm(sample_data(CSS_normalized_phylum_qiime.ps)$Raw_paired_reads ~ sample_data(CSS_normalized_phylum_qiime.ps)$Sample_type))

#
## Ordination testing
#
# We can use analysis of similarities (ANOSIM) to test for significant clustering between samples from 
# different groups. 
# Here, we'll test by the "Group" variable
group_variable = get_variable(CSS_normalized_phylum_qiime.ps,"Sample_type")
# Use the anosim() function
anosim(distance(CSS_normalized_phylum_qiime.ps, "bray"), group_variable)


#
## To try out some other statistical tests, explore the following scripts:
# "scripts/Step4_R_Basic_stats.R", "scripts/Step4_Differential_abundance_testing.R", 
# and for other types of figures, check out "scripts/Step4_R_Basic_stats.R"
#
