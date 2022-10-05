# Load required libraries
##### Set up environment
library("phyloseq")
library("dplyr")
library("ggplot2")
library("data.table")
library("tidyr")
library("forcats")
library("vegan")


# Phyloseq has issue with the phyloseq_to_metagenomeSeq() function. I made the change suggested here:
# https://github.com/joey711/phyloseq/issues/1118

phyloseq_to_metagenomeSeq = function(physeq, ...){
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq)}
  # Coerce count data to vanilla matrix of integers
  countData = round(as(otu_table(physeq), "matrix"), digits=0)
  # Create sample annotation if possible
  if(!is.null(sample_data(physeq,FALSE))){
    ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))  
  } else { 
    ADF = NULL 
  }
  # Create taxa annotation if possible
  if(!is.null(tax_table(physeq,FALSE))){
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq),
                                        physeq@tax_table@.Data,row.names = taxa_names(physeq)))
  } else {
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq),
                                        row.names = taxa_names(physeq)))
  }
  # Create MRexperiment
  if(requireNamespace("metagenomeSeq")){
    mrobj = metagenomeSeq::newMRexperiment(counts = countData, phenoData = ADF, featureData = TDF,...)
    # Calculate normalization factor
    if (sum(colSums(countData > 0) > 1) < ncol(countData)) {
      p = suppressMessages(metagenomeSeq::cumNormStat(mrobj))
    }
    else {
      p = suppressMessages(metagenomeSeq::cumNormStatFast(mrobj))
    }
    mrobj = metagenomeSeq::cumNorm(mrobj, p = p)
    return(mrobj)
  }
}