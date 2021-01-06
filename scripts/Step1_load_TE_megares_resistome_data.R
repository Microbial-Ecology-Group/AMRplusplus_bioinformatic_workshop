#
##
### Loading the resistome results (shotgun reads)
##
#
          
# Sample metadata
sample_metadata <- read.table('./data/megarich_sample_metadata.csv', header=T, sep=',', row.names = 1, quote = "")
sample_metadata # Run the object name to get more information about the file we just loaded
                               
# Load MEGARes counts                                         
TE_amr <- read.table('./data/TE_AMR_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")
# We can convert our amr count object to the otu_table format required for phyloseq
TE_amr <- otu_table(TE_amr, taxa_are_rows = TRUE)
annotations <- read.table('data/megares_full_annotations_v2.0.csv', header=T, row.names=1, sep=",", quote = "")
                                                      
# We can now merge these objects to make a phyloseq object
TE_amr.ps <- merge_phyloseq(TE_amr, tax_table(as.matrix(annotations)), sample_data(sample_metadata))

