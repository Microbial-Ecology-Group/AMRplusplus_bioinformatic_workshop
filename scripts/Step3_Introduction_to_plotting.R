##### Set up environment
# Load required libraries
source('scripts/load_R_packages.R')

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



# We can make a simple boxplot to visualize the values in a column. Y value (Observed) by (~) X value (Sample_type) 
boxplot(sample_metadata$Raw_paired_reads ~ sample_metadata$Sample_type)

#
##
### ggplot layers
##
#

# We can create the initial plot by specifying the data that we'll be using and the x and y variables (columns)
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads))

# A blank ggplot is drawn. Even though the x and y are specified, there are no points or lines in it.
# This is because, ggplot doesn’t assume that you meant a scatterplot or a line chart to be drawn.
# We've only told ggplot what dataset to use and what columns should be used for X and Y axis, but we 
# haven’t explicitly asked it to draw any points.

# Also note that aes() function is used to specify the X and Y axes. 
# That’s because, any information that is part of the source dataframe has to be specified inside the aes() function.

# To create a boxplot, we can add a "geom" layer, geom_boxplot()
# Other options include (and many others)
# geom_point()
# geom_line()
# geom_bar()

# Here's an example of a boxplot
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads)) + 
  geom_boxplot()

# Or we can turn this figure into a graph with points
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads)) + 
  geom_point()

# We can add color based on another variable
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot()


# Specifying the aesthetics "aes()" in the plot vs. in the layers
# Aesthetic mappings can be supplied in the initial ggplot() call, in individual layers, or in some combination of both. 
# All of these calls create the same plot specification:
ggplot(sample_metadata, aes(Sample_type, Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot()
ggplot(sample_metadata, aes(Sample_type, Raw_paired_reads)) + 
  geom_boxplot(aes(color = Sample_type))
ggplot(sample_metadata, aes(Sample_type)) + 
  geom_boxplot(aes(y = Raw_paired_reads, color = Sample_type))
ggplot(sample_metadata) + 
  geom_boxplot(aes(Sample_type, Raw_paired_reads, color = Sample_type))

# Each of the labels has different arguments you can use to change the aesthetics.
ggplot(sample_metadata, aes(Sample_type, Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot(fill= "grey")

# Or we can color based on unique factors
ggplot(sample_metadata, aes(Sample_type, Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot(aes(fill= Sample_type))

# We can also add other layers
ggplot(sample_metadata, aes(Sample_type, Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_point()

# There are many options for layers to add, so be sure to search online for your specific needs
ggplot(sample_metadata, aes(Sample_type, Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1)


#
##
### ggplot labels
##
#

# We can also add layers with labels for the graph using ggtitle(), xlab(), ylab(), and labs
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  ggtitle("Metagenomic sequencing results")

# We can change the x and y axis labels with labs()
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  ggtitle("Metagenomic sequencing results") +
  labs(x = "Treatment Sample_type", y = "raw paired reads")

# Or we can use xlab() and ylab()
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  ggtitle("Metagenomic sequencing results") +
  ylab("raw paired reads") +
  xlab("Treatment Sample_type")

# You can add a caption
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(caption = "(based on data from...)")

# Or a tag, for when labelling a subplot with a letter
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(tag = "A")

#
##
### ggplot themes
##
#

# ggplot has a lot of ways to alter the color schemes and the overall aesthetic of the graph
# we like using "theme()" as an easy way to try other color formatting

# You can alter each specific part of the theme, or use predefined themes such as "theme_classic()" shown below
# Again, there are so many options for what you can modify in the theme(), so look up what you need online and you'll probably find it
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(title = "Metagenomic sequencing results", x = "Treatment Sample_type", y = "raw paired reads") + 
  theme(axis.text.x = element_text( size = 18),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank())

# Here's "theme_classic()" 
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(title = "Metagenomic sequencing results", x = "Treatment Sample_type", y = "raw paired reads") + 
  theme_classic()
 
# Here's "theme_minimal()" 
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(title = "Metagenomic sequencing results", x = "Treatment Sample_type", y = "raw paired reads") + 
  theme_minimal()

# Here's "theme_bw()" 
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(title = "Metagenomic sequencing results", x = "Treatment Sample_type", y = "raw paired reads") + 
  theme_bw()

# Here's "theme_dark()" 
ggplot(sample_metadata, aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(title = "Metagenomic sequencing results", x = "Treatment Sample_type", y = "raw paired reads") + 
  theme_dark()


#
##
### ggplot facets
##
#

# So far, we've only plotted one column at a time from our metadata file.
# With data in the "long" or "melted" form, we get a few more options for plotting more data

# For example, our sample_metadata file has 12 columns.
# We can "gather" these results so that the 4 column with read numbers is made into a unique row for each sample.
# We use the "-" to specify which variables to maintain in the new object, but should not be used for making unique combinations

melted_metadata <- gather(sample_metadata, Raw_paired_reads, Value,
                          -Sample_type , -Sample, -Lot, -Host, -Matrix, -Head, -PREVCAT_A_APLUS, -PREVCAT_ALL,
                          -Sample_block, -shotgun_seq_lane, -shotgun_mean_phred_scores)
melted_phylum_qiime <- psmelt(phylum_qiime.ps)

# Plot the results without ordering the factors
ggplot(melted_phylum_qiime, aes(x = Sample_type , y = Raw_paired_reads, color = Raw_paired_reads)) + 
  geom_boxplot() +
  labs(title = "Metagenomic sequencing results", x = "Sample type", y = "raw paired reads") + 
  theme_classic()

# Note that the factors are ordered alphabetically, but we can also change the order of the factors by using "fct_relevel"
melted_phylum_qiime <- melted_phylum_qiime %>% 
  mutate(Sample_type = fct_relevel(Sample_type,"WWTP", "Swine","Beef","Poultry"))

# Plot the results with ordered factors
ggplot(melted_phylum_qiime, aes(x = Sample_type , y = Raw_paired_reads, color = Raw_paired_reads)) + 
  geom_boxplot() +
  labs(title = "Metagenomic sequencing results", x = "Sample type", y = "raw paired reads") + 
  theme_classic()

#
## Because we are using "melted data", we can also include information about the taxa
# It can also be useful to make "facets" for splitting data into multiple figures
ggplot(melted_phylum_qiime, aes(x = Sample_type , y = Abundance, color = Raw_paired_reads)) + 
  geom_boxplot() +
  labs(title = "Metagenomic sequencing results", x = "Sample type", y = "raw paired reads") + 
  theme_classic() +
  facet_wrap( ~ phylum)

# We can also make sure the scales are free to change depending on the facet.
ggplot(melted_phylum_qiime, aes(x = Sample_type , y = Abundance , color = Raw_paired_reads)) + 
  geom_boxplot() +
  labs(title = "Metagenomic sequencing results", x = "Sample type", y = "raw paired reads") + 
  theme_classic() +
  facet_wrap( ~ phylum, scales = "free")


# Remember, these diversity indices were calculated based on the entire dataset containing all features.
# Often, you'll want to aggregate your counts to other taxonomic levels and calculate diversity indices on those counts.

# Here's an example of how we can do that at the phylum level for the microbiome and drug class level for the resistome


# Aggregate counts to the phylum level, calculate diversity values, and mutate the object to include SeqType and DataType
phylum_kraken.ps <- tax_glom(kraken_microbiome.ps, "phylum")
phylum_kraken_shotgun_diversity_values <- estimate_richness(phylum_kraken.ps, measures = c("Observed","Shannon"))
phylum_kraken_shotgun_diversity_values <- phylum_kraken_shotgun_diversity_values %>%
  mutate(SeqType = "shotgun", DataType = "shotgun microbiome", Sample = row.names(phylum_kraken_shotgun_diversity_values))

# Aggregate counts to the phylum level, calculate diversity values, and mutate the object to include SeqType and DataType
phylum_qiime.ps <- tax_glom(qiime_microbiome.ps, "phylum")
phylum_qiime_diversity_values <- estimate_richness(phylum_qiime.ps, measures = c("Observed","Shannon"))
phylum_qiime_diversity_values <- phylum_qiime_diversity_values %>%
  mutate(SeqType = "16S", DataType = "16S microbiome", Sample = row.names(phylum_qiime_diversity_values))

# Aggregate counts to the drug class level, calculate diversity values, and mutate the object to include SeqType and DataType
class_amr.ps <- tax_glom(amr.ps, "class")
class_amr_shotgun_diversity_values <- estimate_richness(class_amr.ps, measures = c("Observed","Shannon"))
class_amr_shotgun_diversity_values <- class_amr_shotgun_diversity_values %>%
  mutate(SeqType = "shotgun", DataType = "resistome", Sample = row.names(class_amr_shotgun_diversity_values))

# Aggregate counts to the drug class level, calculate diversity values, and mutate the object to include SeqType and DataType
class_TE_amr.ps <- tax_glom(TE_amr.ps, "class")
class_TE_amr_shotgun_diversity_values <- estimate_richness(class_TE_amr.ps, measures = c("Observed","Shannon"))
class_TE_amr_shotgun_diversity_values <- class_TE_amr_shotgun_diversity_values %>%
  mutate(SeqType = "Target enrichment", DataType = "TE resistome", Sample = row.names(class_TE_amr_shotgun_diversity_values))


# Now we can merge these tables based on identical row
phylum_and_class_combined_diversity_values <- bind_rows(phylum_kraken_shotgun_diversity_values, phylum_qiime_diversity_values,class_amr_shotgun_diversity_values,class_TE_amr_shotgun_diversity_values)

# Change the name of the "Sample" column to "Name" so that it matches our other metadata files
phylum_and_class_combined_diversity_values$Name <- phylum_and_class_combined_diversity_values$Sample

# Combine metadata file for MEGARich dataset and the 16S dataset
combined_metadata <- bind_rows(qiime_sample_metadata, sample_metadata)

# Add the metadata file to the new object
phylum_and_class_combined_diversity_values <- left_join(phylum_and_class_combined_diversity_values, combined_metadata, by = "Name")


# Create boxplots for observed features at the phylum and drug class levels
fig7 <- ggplot(phylum_and_class_combined_diversity_values, aes(x = Sample_type, y = Observed, color = Sample_type)) +
  geom_boxplot() +
  labs(title = "Unique phyla and resistome AMR drug classes by data type", x = "Treatment Sample_type", y = "Observed features") + 
  theme_classic() +
  facet_wrap( ~ DataType, scales = "free")
fig7


#
##
### Microbiome plots
##
#

# When plotting the microbiome, we often focus on aggregate counts
# We saw that we can use phyloseq to easily make barplots like this:
plot_bar(phylum_qiime.ps, fill = "phylum")

# The cool thing is that phyloseq plays nicely with ggplot and you use what you learned above to modify these figures
plot_bar(phylum_qiime.ps, fill = "phylum") + 
  facet_wrap(~ Sample_type, scales = "free_x") +
  theme_classic()

# You can also make facets based on the unique phyla in the phyloseq object
# Notice we had to specify that we want "free" scales on both axes
plot_bar(phylum_qiime.ps, fill = "phylum") + 
  facet_wrap(~ phylum, scales = "free") +
  theme_classic()


#
##
### Filtering out low abundance samples
##
#

# As you saw in the figures above, you can end up with too many features that make your figures messy.
# To solve this, we can filter out features that are in low abundance

# only keep OTUs present in greater than 0.5% of all OTUs across all samples.
minTotRelAbun = 0.005
x = taxa_sums(phylum_qiime.ps)
keepTaxa = (x / sum(x)) > minTotRelAbun
length(keepTaxa[keepTaxa ==TRUE])
pruned_phylum_qiime.ps = prune_taxa(keepTaxa, phylum_qiime.ps)
plot_bar(pruned_phylum_qiime.ps, fill = "phylum")


#
##
### Converting data to relative abundances
##
#

# We'll talk more about how to we can account for differences in sequencing depth between samples with count normalization.
# One of the easiest way we can begin to compare the microbiome and resistome composition is by plotting relative abundance
##### Convert OTU abundances to relative abundances
phylum_qiime.ps.rel <- transform_sample_counts(phylum_qiime.ps, function(x) x / sum(x) )

# We can plot these results, notice the y-axis relative abundance plots
plot_bar(phylum_qiime.ps.rel, fill= "phylum")

# We can use the phyloseq object for some of these exploratory figures, but we recommend converting the data into "long" format
phylum_qiime.ps.rel.melt <- psmelt(phylum_qiime.ps.rel)
head(phylum_qiime.ps.rel.melt)

# We can see that we still have 40 unique phyla in this object
unique(phylum_qiime.ps.rel.melt$phylum)

# Lets change the names of phyla to "low abundance phyla" for features present at a relative
# proportion less than 0.05 across all samples
phylum_qiime.ps.rel.melt <- phylum_qiime.ps.rel.melt %>%
  group_by(phylum) %>%
  mutate(mean_phylum_rel_abundance = mean(Abundance))

phylum_qiime.ps.rel.melt$phylum <- as.character(phylum_qiime.ps.rel.melt$phylum)
phylum_qiime.ps.rel.melt$mean_phylum_rel_abundance <- as.numeric(phylum_qiime.ps.rel.melt$mean_phylum_rel_abundance)

phylum_qiime.ps.rel.melt$phylum[phylum_qiime.ps.rel.melt$mean_phylum_rel_abundance < 0.005] <- "Low abundance phyla"

##### Plot phyla relative abundances
ggplot(phylum_qiime.ps.rel.melt, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Sample_type, scales = "free") +
  theme_classic()


## We can use the data in long form with the dplyr functions we learned to summarize the results

# Mean relative abundance of phyla across all samples
phylum_qiime.ps.rel.melt %>%
  group_by(phylum) %>%
  summarize(mean_rel_abundance = mean(Abundance)) %>%
  arrange(-mean_rel_abundance)

# Mean relative abundance of phyla by treatment Sample_type
phylum_qiime.ps.rel.melt %>%
  group_by(Sample_type, phylum) %>%
  summarize(mean_rel_abundance = mean(Abundance)) %>%
  arrange(-mean_rel_abundance)



#
##
### Rarefaction plots
##
#

# Rarefaction is a technique used to estimate how many features (labeled "species") would be identified
# at different sequencing depths (labeled "Sample Size"). 

# Generally, if the curves are relatively flat at the actual sample size, this suggests adequate sequencing depth.
# Alternatively, steep curves suggest that increasing the sequencing depth 

# Rarefaction for the qiime2 microbiome at the phylum level
# Notice, that we had to get the "otu_table()" from the phyloseq object, then we had use "t()" to transpose the table
phylum_qiime_rarefaction <- rarecurve(t(otu_table(phylum_qiime.ps)), step = 200, se = TRUE)
# Let's see if the rarefaction looks different at the genus level
genus_qiime.ps <- tax_glom(qiime_microbiome.ps, "genus")
genus_qiime_rarefaction <- rarecurve(t(otu_table(genus_qiime.ps)), step = 200, se = TRUE)

# Rarefaction for the kraken microbiome at the phylum level
phylum_kraken_rarefaction <- rarecurve(t(otu_table(phylum_kraken.ps)), step = 200, se = TRUE)
# Let's see if the rarefaction looks different at the genus level
genus_kraken.ps <- tax_glom(kraken_microbiome.ps, "genus")
genus_kraken_rarefaction <- rarecurve(t(otu_table(genus_kraken.ps)), step = 200, se = TRUE)

# Rarefaction for the resistome at the AMR drug class level
class_amr_rarefaction <- rarecurve(t(otu_table(class_amr.ps)), step = 200, se = TRUE)
# Let's see if the rarefaction looks different at the AMR gene group level
group_amr.ps <- tax_glom(amr.ps, "group")
group_amr_rarefaction <- rarecurve(t(otu_table(group_amr.ps)), step = 200, se = TRUE)




