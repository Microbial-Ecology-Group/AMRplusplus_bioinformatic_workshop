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
### Load qiime2 microbiome results (16S reads)
##
#

# We import the .biom as in Lesson 1 step 4
qiime_microbiome <- import_biom("data/Exported_16S_qiime2_results/16S_asv_table.biom")
# Load the taxonomy file from qiime2
taxa <- read.table("data/Exported_16S_qiime2_results/taxonomy.tsv", header=T, row.names=1, sep='\t', quote = "")
row.names(taxa) <- paste(row.names(taxa),taxa[,1], sep= '; ')
taxa.dt <- data.table(id=rownames(taxa)) # we'll make a column with the name "id"
taxa.dt[, c('feature',
            'kingdom',
            'phylum',
            'class',
            'order',
            'family',
            'genus',
            'species') := tstrsplit(id, '; ', type.convert = TRUE, fixed = TRUE)]

taxa.df <- as.data.frame(taxa.dt)
taxa.df <- within(taxa.df, rm(id))
row.names(taxa.df) <- taxa.df$feature
taxa.df <- within(taxa.df, rm(feature))

qiime_microbiome_phylo_tree <- read_tree("./data/Exported_16S_qiime2_results/tree.nwk")

qiime_microbiome.ps <- merge_phyloseq(qiime_microbiome, phy_tree(qiime_microbiome_phylo_tree), tax_table(as.matrix(taxa.df)), sample_data(sample_metadata))

# Estimating richness and diversity using the easy-to-use function estimate_richness()
qiime_microbiome_16S_diversity_values <- estimate_richness(qiime_microbiome.ps)

