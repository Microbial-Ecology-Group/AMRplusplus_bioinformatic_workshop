##### Set up environment
library("phyloseq")
library("dplyr")
library("ggplot2")
library("data.table")
library("tidyr")
library("forcats")
library("vegan")


#
##
### Loading the resistome results (shotgun reads)
##
#
                                         
# Load MEGARes counts                                         
amr <- read.table('./data/shotgun_AMR_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")
# We can convert our amr count object to the otu_table format required for phyloseq
amr <- otu_table(amr, taxa_are_rows = TRUE)
annotations <- read.table('data/megares_full_annotations_v2.0.csv', header=T, row.names=1, sep=",", quote = "")
                                                      
# We can now merge these objects to make a phyloseq object
amr.ps <- merge_phyloseq(amr, tax_table(as.matrix(annotations)), sample_data(sample_metadata))

# Estimating richness and diversity using the easy-to-use function estimate_richness()
amr_shotgun_diversity_values <- estimate_richness(amr.ps)


