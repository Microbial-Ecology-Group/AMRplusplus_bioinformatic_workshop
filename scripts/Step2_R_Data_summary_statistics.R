##### Set up environment
library("phyloseq")
library("dplyr")
library("ggplot2")
library("data.table")

# Load qiime2 data
source("scripts/Step1_load_qiime2_microbiome_data.R")
# Check out this new phyloseq object
qiime_microbiome.ps # This phyloseq object contains the microbiome results, taxa table, and metadata


# We can access the individual components like this:
# Look at this site for a few more examples: https://joey711.github.io/phyloseq/import-data.html
otu_table(qiime_microbiome.ps)
tax_table(qiime_microbiome.ps)
sample_data(qiime_microbiome.ps)
phy_tree(qiime_microbiome.ps)

# Here are some other functions that provide information about the phyloseq object
# Just the sample names
sample_names(qiime_microbiome.ps)
# Taxonomic rank names
rank_names(qiime_microbiome.ps)
# Variables in the metadata
sample_variables(qiime_microbiome.ps)
# Easy sums of counts by sample
sample_sums(qiime_microbiome.ps)

# How can we get the total counts by adding across all samples?





#
##
### Metadata exploration
##
#
#


# Next, we might want to make multiple calculations on our metadata variables
# Here, we can use the "dplyr" package with the "%>%" function and the "summarize()" function as shown below:
sample_data(qiime_microbiome.ps) %>%
  summarize(total_16S_counts = sum(Raw_paired_reads), mean_16S_counts = mean(Raw_paired_reads),
            min_16S_counts = min(Raw_paired_reads),max_16S_counts = max(Raw_paired_reads),
            median_16S_counts = median(Raw_paired_reads))

# We can group our summary statistics by specific metadata variables by adding the "group_by()" function
# Ignore the warning about the R class, this shouldn't affect your downstream analysis
sample_data(qiime_microbiome.ps) %>%
  group_by(Sample_type) %>% 
  summarize(total_16S_counts = sum(Raw_paired_reads), mean_16S_counts = mean(Raw_paired_reads),
            min_16S_counts = min(Raw_paired_reads),max_16S_counts = max(Raw_paired_reads),
            median_16S_counts = median(Raw_paired_reads))



#
## Agglomerate ASV counts to different taxonomic levels
#

# Using tax_glom(), we can easily aggregate counts to different levels (taxonomic or AMR annotation levels)
qiime_phylum.ps <- tax_glom(qiime_microbiome.ps, "phylum")
qiime_phylum.ps

# We can get sample counts at the phylum level.
sum(sample_sums(qiime_phylum.ps))

# By summing counts at each taxonomic level, we can observe the reduction in sample counts in lower taxonomic levels.
# Note the stop sign on the right hand side of the Console panel, when you see this it means R is running a command.
# You can click on the stop sign to stop the run.
species.ps <- tax_glom(qiime_microbiome.ps, "species") # Depending on your computer and file size, this can take a couple of seconds to a few minutes

# We can get sample counts at the species level.
sum(sample_sums(species.ps))

# We can calculate the percentage of counts at the species level out of classified at the phylum level.
sum(sample_sums(species.ps))  /  sum(sample_sums(qiime_phylum.ps)) * 100

# In this data, we had a high percentage of ASVs classified down to the species level. 
# This can vary widely by dataset and is an important consideration.
# For example, we might not emphasize results from the species-level counts in a manuscript if they only made up ~ 50% of counts.



#
## Visualize the microbiome composition
#

# We can observe the counts using a simple bar chart using plot_bar()
# We have to specify that we'll fill the bar charts using the phylum counts ()
plot_bar(qiime_microbiome.ps)


plot_bar(qiime_phylum.ps, fill = "phylum")

# We can specify a "facet_grid" to organize the plot
plot_bar(qiime_phylum.ps, fill = "phylum", facet_grid = "Sample_type")


# We group counts by the treatement Sample_type
plot_bar(qiime_phylum.ps, fill = "phylum", x = "Sample_type")


##### Convert OTU abundances to relative abundances
qiime_phylum.ps.rel <- transform_sample_counts(qiime_phylum.ps, function(x) x / sum(x) )

# We can plot these results, notice the y-axis relative abundance plots
plot_bar(qiime_phylum.ps.rel, fill= "phylum")





#
## Microbiome diversity indices 
#



# Alpha-diversity indices are a common measurement to summarize the composition of the microbiome and resistome.
# "Richness" or "Observed" simply describes the number of unique taxa identified in a sample
# On the other hand, we use "Shannon's index" or "Inverse Simpsons's index" to describe "evenness", or
# distribution of counts among taxa in a sample

# Find more infomation about diversity indices here: http://www.jmb.or.kr/submission/Journal/027/JMB027-12-02_FDOC_2.pdf

# Estimating richness and diversity using the easy-to-use function estimate_richness()
qiime_microbiome_16S_diversity_values
qiime_microbiome_16S_diversity_values$Name <- row.names(qiime_microbiome_16S_diversity_values)

# Remember, that with phyloseq objects, you have to convert to matrix first, then to dataframe
qiime_sample_metadata <- as.data.frame(as(sample_data(qiime_microbiome.ps), "matrix"))

# Now, we use left_join() to add the sample_metadata to the combined_diversity_values object
qiime_microbiome_16S_diversity_values <- left_join(qiime_sample_metadata,qiime_microbiome_16S_diversity_values, by = "Name")

# Need to provide sample_names() to the newly formed metadata object
sample_names(qiime_microbiome_16S_diversity_values) <- qiime_microbiome_16S_diversity_values$Name
# Merge updated metadata values back into phyloseq object
qiime_microbiome.ps <- merge_phyloseq(qiime_microbiome, phy_tree(qiime_microbiome_phylo_tree), tax_table(as.matrix(taxa.df)), sample_data(sample_metadata))


# We can do the same for the phylum counts
qiime_phylum_diversity <- estimate_richness(qiime_phylum.ps) # Notice the error you get.

# Microbiome papers usually report diversity indices at the ASV level (before aggregation), but also for different taxa levels

# We can easily plot these values using the plot_richness() function, we'll just pick 3 commonly used alpha diversity indices
plot_richness(qiime_phylum.ps, x = "Sample_type", color = "Sample_type", measures = c("Observed", "Shannon", "InvSimpson")) 

# Try again for just a few diversity measures
qiime_phylum_diversity <- estimate_richness(qiime_phylum.ps, measures = c("Observed", "Shannon", "InvSimpson"))
qiime_phylum_diversity

# We can then modify this figure using ggplot2 functions like geom_boxplot() by using the "+" sign 
# We'll explore further options in later steps
plot_richness(qiime_phylum.ps,x = "Sample_type", color = "Sample_type", measures = c("Observed", "Shannon", "InvSimpson")) + 
  geom_boxplot()





#
##
### Export phyloseq object to "melted" format for summary statistics
##
#

# While phyloseq has built-in functions for easy plotting, we might want to export these counts for easier use with ggplot2

##### Convert phyloseq object with relative abundance to data.frame
qiime_phylum.ps.rel.melt <- psmelt(qiime_phylum.ps.rel)

head(qiime_phylum.ps.rel.melt)

##### Plot phyla abundances
plot.phylum.rel.bar <- ggplot(qiime_phylum.ps.rel.melt, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Sample_type, scales = "free") +
  theme_bw()

plot.phylum.rel.bar


#
##
### Use this area below to try the methods you just learned, but using another 
##
#