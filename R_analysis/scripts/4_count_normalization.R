# load previous work
source("Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/scripts/3_alpha_diversity.R")


#### NORMALIZATION WITH CUMULATIVE SUM SCALING (CSS) ####

# CSS is currently our preferred count normalization method
# CSS calculates the quantile of the count distribution of samples where they all should be roughly 
# equivalent and independent of each other up to this quantile under the assumption that, 
# at this range, counts are derived from a common distribution.
# CSS normalization is implemented in the "metagenomeSeq" R library

# First, we convert the phyloseq object to metagenomeSeq
microbiome.metaseq <- phyloseq_to_metagenomeSeq(microbiome.ps) 
amr.metaseq <- phyloseq_to_metagenomeSeq(amr.ps)

# let's look at them
microbiome.metaseq
amr.metaseq

# Perform the actual normalization
microbiome_css.metaseq <- cumNorm(microbiome.metaseq)
amr_css.metaseq <- cumNorm(amr.metaseq)

# Like phyloseq, metagenomeSeq has it's own functions for accessing their data.
# Here, we need to use MRcounts() and re-make the phyloseq object with the normalized counts
CSS_microbiome_counts <- MRcounts(microbiome.metaseq, norm = TRUE)
CSS_amr_counts <- MRcounts(amr.metaseq, norm = TRUE)

# Use the new counts and merge with components from our original phyloseq object.
microbiome_css.ps <- merge_phyloseq(otu_table(CSS_microbiome_counts, taxa_are_rows = TRUE),sample_data(microbiome.ps),tax_table(microbiome.ps), phy_tree(microbiome.ps))
amr_css.ps <- merge_phyloseq(otu_table(CSS_amr_counts, taxa_are_rows = TRUE),sample_data(amr.ps),tax_table(amr.ps))

# aggregate the normalized counts to the phylum and class level
microbiome_ccs_phylum.ps <- tax_glom(microbiome_css.ps, taxrank = "Phylum", NArm = FALSE)
amr_css_class.ps <- tax_glom(amr_css.ps, taxrank = "class")

#### LET'S PLOT THE PHYLA BY SAMPLE TYPE
plot_microbiome_css_phylum <- plot_bar(microbiome_ccs_phylum.ps, fill = "Phylum") +
  facet_wrap(~Sample_type, scales = "free_x")
plot_microbiome_css_phylum

plot_amr_css_class <- plot_bar(amr_css_class.ps, fill = "class") +
  facet_wrap(~Sample_type, scales = "free_x")
plot_amr_css_class
# Notice the y-axis values are counts and not proportions.
