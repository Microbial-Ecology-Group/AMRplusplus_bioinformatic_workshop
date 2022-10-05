# Edit full dataset
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
### Load in the metadata file
##
#


sample_metadata <- read.table('./data/proj7_sample_metadata.csv', header=T, sep=',', row.names = 1, quote = "")
sample_metadata # Run the object name to get more information about the file we just loaded

#
##
### Loading the kraken2 microbiome results
##
#

# Load the kraken count table                                        
kraken_microbiome <- read.table('./data/full_kraken_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")

# Convert to format that phyloseq likes with otu_table()                                      
kraken_microbiome <- otu_table(kraken_microbiome, taxa_are_rows = TRUE)

# Repeat similar steps to what we did with the qiime2 taxonomy
kraken_taxonomy <- data.table(id=rownames(kraken_microbiome))
kraken_taxonomy[, c('domain',
                    'kingdom',
                    'phylum',
                    'class',
                    'order',
                    'family',
                    'genus',
                    'species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]

# Conver to data.frame
kraken_taxonomy <- as.data.frame(kraken_taxonomy)
# Use the id variable to rename the row.names
row.names(kraken_taxonomy) <- kraken_taxonomy$id
# Remove the "id" column
kraken_taxonomy <- within(kraken_taxonomy, rm(id))

# Create kraken phyloseq object
kraken_microbiome.ps <- merge_phyloseq(kraken_microbiome, tax_table(as.matrix(kraken_taxonomy)), sample_data(sample_metadata))

# Subset for just shotgun samples
kraken_microbiome.ps <- subset_samples(kraken_microbiome.ps, Library_Type=="shotgun")

# Print out count matrix
write.csv(otu_table(kraken_microbiome),"data/shotgun_kraken_analytic_matrix.csv", quote = FALSE)


#
##
### Loading MEGARes resistome results
##
#

# Load MEGARes counts                                         
amr <- read.table('./data/full_AMR_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")
# We can convert our amr count object to the otu_table format required for phyloseq
amr <- otu_table(amr, taxa_are_rows = TRUE)
annotations <- read.table('data/megares_full_annotations_v2.0.csv', header=T, row.names=1, sep=",", quote = "")

# We can now merge these objects to make a phyloseq object
amr.ps <- merge_phyloseq(amr, tax_table(as.matrix(annotations)), sample_data(sample_metadata))

# Estimating richness and diversity using the easy-to-use function estimate_richness()
amr_shotgun_diversity_values <- estimate_richness(amr.ps)

# Subset for just shotgun samples
shotgun_amr.ps <- subset_samples(amr.ps, Library_Type=="shotgun")
TE_amr.ps <- subset_samples(amr.ps, Library_Type=="TE")


# Print out count matrix
write.csv(otu_table(shotgun_amr.ps),"data/shotgun_AMR_analytic_matrix.csv", quote = FALSE)
write.csv(otu_table(TE_amr.ps),"data/TE_AMR_analytic_matrix.csv", quote = FALSE)



