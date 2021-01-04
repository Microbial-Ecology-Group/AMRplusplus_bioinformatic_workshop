#
##
### Loading the kraken2 microbiome results (shotgun reads)
##
#

# Load the kraken count table                                        
kraken_microbiome <- read.table('./data/shotgun_kraken_analytic_matrix.csv', header=T, row.names=1, sep=',')

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
# Remove the "id" and "kingdom" columns
kraken_taxonomy <- within(kraken_taxonomy, rm(id))
kraken_taxonomy <- within(kraken_taxonomy, rm(kingdom))

# Create kraken phyloseq object
kraken_microbiome.ps <- merge_phyloseq(kraken_microbiome, tax_table(as.matrix(kraken_taxonomy)),sample_data(sample_metadata))

# Estimating richness and diversity using the easy-to-use function estimate_richness()
microbiome_shotgun_diversity_values <- estimate_richness(kraken_microbiome.ps)

