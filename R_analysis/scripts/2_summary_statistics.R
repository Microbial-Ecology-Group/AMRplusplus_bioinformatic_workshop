#load R packages
source("Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/scripts/load_packages.R")
# load TE and 16S data
source("Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/scripts/1_load_data.R")

# look at the phyloseq objects
amr.ps # this contains the TE-AMR counts, taxonomy, and metadata
microbiome.ps # this one the 16S-ASV counts, taxonomy, tree, rep-seqs, and metadata


# We can access the individual components like this:
# Look at this site for a few more examples: https://joey711.github.io/phyloseq/import-data.html
# counts
otu_table(microbiome.ps)
otu_table(amr.ps)
# taxonomy
tax_table(microbiome.ps)
tax_table(amr.ps)
# metadata
sample_data(microbiome.ps)
sample_data(amr.ps)
# tree (16S only)
phy_tree(microbiome.ps)

# Here are some other functions that provide information about the phyloseq object
# Just the sample names
sample_names(microbiome.ps)
sample_names(amr.ps)
# Taxonomic rank names
rank_names(microbiome.ps)
rank_names(amr.ps)
# Variables in the metadata
sample_variables(microbiome.ps)
sample_variables(amr.ps)

# Easy sums of counts by sample
sample_sums(microbiome.ps)
sample_sums(amr.ps)

# How can we get the total *ASVs* by adding across all samples?
sum(sample_sums(microbiome.ps))
mean(sample_sums(microbiome.ps))
min(sample_sums(microbiome.ps))

sum(sample_sums(amr.ps))
mean(sample_sums(amr.ps))
min(sample_sums(amr.ps))


#### METADATA EXPLORATION ####

# Next, we might want to make multiple calculations on our metadata variables
# Here, we can use the "dplyr" package with the "%>%" function and the "summarize()" function as shown below:
# Calculate the numbers of *reads* 
sample_data(microbiome.ps) %>%
  summarize(total_16S_counts = sum(Raw_paired_reads), mean_16S_counts = mean(Raw_paired_reads),
            min_16S_counts = min(Raw_paired_reads),max_16S_counts = max(Raw_paired_reads),
            median_16S_counts = median(Raw_paired_reads))

sample_data(amr.ps) %>%
  summarize(total_ARG_counts = sum(Raw_paired_reads), mean_ARG_counts = mean(Raw_paired_reads),
            min_ARG_counts = min(Raw_paired_reads),max_ARG_counts = max(Raw_paired_reads),
            median_ARG_counts = median(Raw_paired_reads))

# We can group our summary statistics by specific metadata variables by adding the "group_by()" function
sample_data(microbiome.ps) %>%
  group_by(Sample_type) %>% 
  summarize(total_16S_counts = sum(Raw_paired_reads), mean_16S_counts = mean(Raw_paired_reads),
            min_16S_counts = min(Raw_paired_reads),max_16S_counts = max(Raw_paired_reads),
            median_16S_counts = median(Raw_paired_reads))

sample_data(amr.ps) %>%
  group_by(Sample_type) %>% 
  summarize(total_ARG_counts = sum(Raw_paired_reads), mean_ARG_counts = mean(Raw_paired_reads),
            min_ARG_counts = min(Raw_paired_reads),max_ARG_counts = max(Raw_paired_reads),
            median_ARG_counts = median(Raw_paired_reads))

#### plotting and comparing sampling depth ####
## 16S
sequencing_depth_16S_boxplot = ggplot(sample_data(microbiome.ps), aes(x= Sample_type, y= Raw_paired_reads, fill = Sample_type)) +
  theme_bw() +
  geom_boxplot() +
  geom_point()
sequencing_depth_16S_boxplot

pairwise.wilcox.test(sample_data(microbiome.ps)$Raw_paired_reads, sample_data(microbiome.ps)$Sample_type) #NS

## TE
sequencing_depth_amr_boxplot <- ggplot(sample_data(amr.ps), aes(x= Sample_type, y= Raw_paired_reads, fill = Sample_type)) +
  theme_bw() +
  geom_boxplot() +
  geom_point()
sequencing_depth_amr_boxplot

# test it
pairwise.wilcox.test(sample_data(microbiome.ps)$Raw_paired_reads, sample_data(microbiome.ps)$Sample_type) #NS

pairwise.wilcox.test(sample_data(amr.ps)$Raw_paired_reads, sample_data(microbiome.ps)$Sample_type) #NS


#### AGGLOMERATE TO DIFFERENT TAXONOMIC LEVELS ####

# Using tax_glom(), we can easily aggregate counts to different levels (taxonomic or AMR annotation levels)
microbiome_phylum.ps <- tax_glom(microbiome.ps, taxrank = "Phylum", NArm = T)
amr_class.ps <- tax_glom(amr.ps, taxrank = "class")

microbiome_phylum.ps # 37 phyla across all samples
amr_class.ps # 48 classes across all samples
