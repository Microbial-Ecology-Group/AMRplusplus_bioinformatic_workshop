# load packages
source("Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/scripts/load_packages.R")

######## LOADING DATA ########

#### TE DATA ####
# read in AMR count matrix
amr_data <- read.table('Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/data/TE_AMR_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")
# convert this into an 'otu_table' format required for phyloseq
amr_data <- otu_table(amr_data, taxa_are_rows = T)

#read in gene annotations
annotations <- read.table('Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/data/megares_modified_annotations_v2.00.csv', header=T, row.names=1, sep=",", quote = "")
#convert this into a 'taxonomy table' format required for phyloseq
annotations <- tax_table(as.matrix(annotations))

# read in TE metadata
amr_metadata <- read.table('Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/data/TE_metadata.csv', header=T, sep=',', row.names = 1, quote = "")
# convert to 'sample_data' format for phyloseq
amr_metadata <- sample_data(amr_metadata)

# merge the annotations, the count matrix, and metadata into a phyloseq object
amr.ps <- merge_phyloseq(amr_data, annotations, amr_metadata)


#### 16S DATA ####
# read in ASV table, tree, representative sequences
microbiome_data <- import_biom('Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/data/table-with-taxonomy.biom','Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/data/tree.nwk','Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/data/dna-sequences.fasta')


# read in metadata
microbiome_metadata <- read.table('Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/data/16S_metadata.csv', header=T, sep=',', row.names = 1, quote = "")
# convert to 'sample_data' format for phyloseq
microbiome_metadata <- sample_data(microbiome_metadata)

# merge the ASVs, tree, rep-seqs, and metadata into phyloseq object
microbiome.ps <- merge_phyloseq(microbiome_data, microbiome_metadata)
# convert columns of taxonomy table to match taxonomic levels
colnames(tax_table(microbiome.ps)) <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus")

# changing the SILVA style naming (k__Bacteria, etc.)
tax.data <- data.frame(tax_table(microbiome.ps)) # extract the taxonomy table as a data frame
tax.data.names <- changeSILVAtaxa(tax.data) # this gets rid of the SILVA style
# now to change the NAs to a better naming scheme
for (i in 1:7){ tax.data.names[,i] <- as.character(tax.data.names[,i])} # converting all columns to characters
tax.data.names[is.na(tax.data.names)] <- "" # replacing the NAs with an empty string
# now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == ""){
    domain <- paste("unclassified ", tax.data.names[i,1], sep = "")
    tax.data.names[i, 2:7] <- domain
  } else if (tax.data.names[i,3] == ""){
    kingdom <- paste("unclassified ", tax.data.names[i,2], sep = "")
    tax.data.names[i, 3:7] <- kingdom
  } else if (tax.data.names[i,4] == ""){
    phylum <- paste("unclassified ", tax.data.names[i,3], sep = "")
    tax.data.names[i, 4:7] <- phylum
  } else if (tax.data.names[i,5] == ""){
    class <- paste("unclassified ", tax.data.names[i,4], sep = "")
    tax.data.names[i, 5:7] <- class
  } else if (tax.data.names[i,6] == ""){
    order <- paste("unclassified ", tax.data.names[i,5], sep = "")
    tax.data.names[i, 6:7] <- order
  } else if (tax.data.names[i,7] == ""){
    tax.data.names$Genus[i] <- paste("unclassified ",tax.data.names$Family[i], sep = "")
  }
}
tax_table(microbiome.ps) <- as.matrix(tax.data.names) # re-insert the taxonomy table into the phyloseq object
tax_table(microbiome.ps)