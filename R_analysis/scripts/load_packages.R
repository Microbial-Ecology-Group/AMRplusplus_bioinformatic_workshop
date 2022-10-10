#### LOAD LIBRARIES ####
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(vegan)
library(data.table)
library(metagenomeSeq)
library(stringr)
library(ggdendro)
library(cowplot)
library(pairwiseAdonis)


#### LOAD FUNCTION TO CHANGE SILVA FORMATTED TAXA NAMES####
changeSILVAtaxa <- function(x) {
  # remove the D__ etc...
  tax.clean <- data.frame(row.names = row.names(x),
                          Domain = str_replace(x[,1], "d__",""),
                          Kingdom = str_replace(x[,2], "k__",""),
                          Phylum = str_replace(x[,3], "p__",""),
                          Class = str_replace(x[,4], "c__",""),
                          Order = str_replace(x[,5], "o__",""),
                          Family = str_replace(x[,6], "f__",""),
                          Genus = str_replace(x[,7], "g__",""),
                          stringsAsFactors = FALSE)
}
