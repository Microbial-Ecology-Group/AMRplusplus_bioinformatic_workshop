# Load required libraries
source('scripts/load_R_packages.R')
#
##
### Load in the metadata file
##
#


sample_metadata <- read.table('./data/megarich_16S_sample_metadata.csv', header=T, sep=',', row.names = 1, quote = "")
sample_metadata # Run the object name to get more information about the file we just loaded

# Change name of the first column
#colnames(sample_metadata)[1] <- "Sample_name"



#
##
### Load qiime2 microbiome results (16S reads)
##
#

# We import the .biom as in Lesson 1 step 4
qiime_microbiome <- import_biom("data/Megarich_exported_16S_qiime2_results/otu_table_json.biom")
# Load the taxonomy file from qiime2
taxa <- read.table("data/Megarich_exported_16S_qiime2_results/taxonomy.tsv", header=T, row.names=1, sep='\t', quote = "")
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

# This could be done in a simpler way, but we need to clean up the taxa.df before using it with phyloseq
taxa.df <- as.data.frame(taxa.dt)
taxa.df <- within(taxa.df, rm(id))
row.names(taxa.df) <- taxa.df$feature
taxa.df <- within(taxa.df, rm(feature))

qiime_microbiome_phylo_tree <- read_tree("./data/Megarich_exported_16S_qiime2_results/tree.nwk")

qiime_microbiome.ps <- merge_phyloseq(qiime_microbiome, phy_tree(qiime_microbiome_phylo_tree), tax_table(as.matrix(taxa.df)), sample_data(sample_metadata))
