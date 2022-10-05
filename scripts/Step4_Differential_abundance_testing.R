# Load required libraries
source('scripts/load_R_packages.R')
# Load kraken
source("scripts/Step1_load_kraken_microbiome_data.R") 

# Filtering low abundance phyla prior to normalization.
# For this example, we'll just choose an arbitrary value of 2 as the threshold for
# for filtering out ASVs with low abundance.
filtered_microbiome.ps = filter_taxa(kraken_microbiome.ps, function(x) sum(x) > 2, TRUE)

# We'll use CSS normalized values, so we need to convert to a metagenomeSeq object
filtered_microbiome.metaseq <- phyloseq_to_metagenomeSeq(filtered_microbiome.ps)

# This is where we perform the CSS normalization
microbiome.metaseq <- cumNorm(filtered_microbiome.metaseq)

# Like phyloseq, metagenomeSeq has it's own functions for accessing their data.
# Here, we need to use MRcounts() and re-make the phyloseq object with the normalized counts
CSS_microbiome_counts <- MRcounts(microbiome.metaseq, norm = TRUE)
# Use the new counts and merge with components from our original phyloseq object.
CSS_normalized_kraken.ps <- merge_phyloseq(otu_table(CSS_microbiome_counts, taxa_are_rows = TRUE),sample_data(kraken_microbiome.ps),tax_table(kraken_microbiome.ps))

# Aggregate counts to phylum
CSS_normalized_phylum_kraken.ps <- tax_glom(CSS_normalized_kraken.ps, "phylum")

# For the functions below to work, we'll need to round the count values
otu_table(CSS_normalized_phylum_kraken.ps) <- round(otu_table(CSS_normalized_phylum_kraken.ps))

# As we've seen in the previous scripts, you would now use these normalized counts for most subsequent analysis
# In this script, we'll need to convert our phylum counts back to an MRexperiment (metagenomeSeq package)
# to use the zero-inflated Gaussian model (ZIG) for differential abundance testing.

# We'll use the package below which provides useful functions for differential abundance testing
library(devtools)
install_github("Russel88/DAtest")
library(DAtest)

# We can use the testDA() function to try out multiple differential abundance tests.
# Check out this website for all other options: https://labs.cd2h.org/gitforager/repository/repository.jsp?id=95257642

# We can use the handy "preDA()" function to add more filtering steps.
# Here, lets just clean up our data and remove any features that aren't in at least 1 sample.
# not needed in this case
CSS_normalized_phylum_kraken.ps <- preDA(CSS_normalized_phylum_kraken.ps, min.samples = 1)

# Using the "testDA" function, we can create the model and test a few different
# methods simultaneously.
# Here, we'll just test the ZIG model and a Wilcoxon test.

DA_results <- testDA(CSS_normalized_phylum_kraken.ps, tests = c("zig", "wil") ,
                     predictor = "Sample_type")

# You can see that the result is actually a list of 20 results from permutations for each test
DA_results
str(DA_results)
# We can view comparisons between the two tests
plot(DA_results)
# Plot a histogram of p-values
plot(DA_results,p = TRUE)

# We can see the results of each permutation by accessing items in the list
DA_results$results[[1]]
# Here's the zig model for the first permutation
DA_results$results[[1]]$zig
# Wilcoxon results for the first permutation
DA_results$results[[1]]$wil

# For our example, we'll stick with using the ZIG model which is implemented in
# the "DA.zig()" function of the "DAtest" library we just installed.
# Please note that all of the actual functions associated with the ZIG model
# are actually a part of the "metagenomeSeq" library
DAzig_model <- DA.zig(CSS_normalized_phylum_kraken.ps, predictor = "Sample_type",
                    allResults = TRUE, p.adj = "bonferroni")
# View results from ZIG model
MRfulltable(DAzig_model)

# Plot association between abundance of a feature and predictor, modified if paired and covars are available
featurePlot(CSS_normalized_phylum_kraken.ps, predictor = "Sample_type",feature = "Bacteria|NA|Actinobacteria|Actinobacteria|Corynebacteriales|Corynebacteriaceae|Corynebacterium|Corynebacterium stationis")


#
##
### Using metagenomeSeq
##
#

# We debated showing you this next code because it adds a lot of other functions you aren't 
# familiar with, but we wanted to make sure you had a template to use in the future.

# Above, we used the "testDA" library because it's functions wraps a bunch of other useful
# functions and allows us to compare the results from different tests.

# We then used the "DA.zig()" function to create the ZIG model, but this is also just a "wrapper"
# for the functions in metagenomeSeq and can limit your control over the model. 

# Below we'll show you the native functions used for creating a ZIG model.

# We can use metagenomeSeq's function aggTax to aggregate counts as in phyloseq
phylum_microbiome.metaseq <- aggTax(microbiome.metaseq, lvl = "phylum")

# metagenomeSeq also has functions for filtering data
filtered_phylum_microbiome.metaseq <- filterData(phylum_microbiome.metaseq, present = 2)

# Make the "zero" model with library size of the raw data
zero_mod <- model.matrix(~0+log(libSize(microbiome.metaseq)))
# Make model with "Sample_type" variable
Sample_type <- pData(phylum_microbiome.metaseq)$Sample_type
design_group = model.matrix(~0 + Sample_type)

# We still need to use cumNorm even thought we aren't using the normalization factor because of the "useCSSoffset = FALSE)                                     
filtered_phylum_microbiome.metaseq <- cumNorm(filtered_phylum_microbiome.metaseq)
                                     
# Create ZIG model                                     
zig_model <- fitZig(obj= filtered_phylum_microbiome.metaseq, mod = design_group, zeroMod=zero_mod, useCSSoffset = FALSE)

# Use Ebayes to adjust model fit
ebayes_zig_model <- eBayes(zig_model@fit)
ebayes_zig_model

zigFit_Group = zig_model@fit
finalMod_Group = ebayes_zig_model$design

contrast_Group = makeContrasts(Sample_typeBeef - Sample_typePoultry,Sample_typeBeef - Sample_typeSwine,
                               Sample_typeBeef - Sample_typeWWTP, Sample_typePoultry - Sample_typeSwine,
                               Sample_typePoultry - Sample_typeWWTP, Sample_typeSwine - Sample_typeWWTP,
                               levels=finalMod_Group)


zigFit_contrasts = contrasts.fit(zigFit_Group, contrast_Group)
EB_zigFit_contrasts = eBayes(zigFit_contrasts)

# make table of results, you can output these results with "write.csv()"
# Explore this table for the results.
full_table_EB_zigFit_contrasts <- topTable(EB_zigFit_contrasts, adjust.method="BH",number = 1000)
# This table has all of the results from the contrasts made above. 
full_table_EB_zigFit_contrasts

# To pick only a comparison between certain coefficient pairs, use the "coef" flag and specify the pair you want to use based on it's position in the table (1,2,3, etc)
table_EB_zigFit_contrasts <- topTable(EB_zigFit_contrasts, adjust.method="BH",number = 1000, coef = 1)
table_EB_zigFit_contrasts

#
## Make simple volcano plots of ZIG results
#
# We'll go over other better ways to plot your ZIG results in Lesson 3, step 5.

# Create sizes for the points on the plot based on "AveExpr"
radius_phylum_node <- sqrt(table_EB_zigFit_contrasts$AveExpr/pi)
# Create the volcano plot, circle size by AveExpr and color by p-value
sign.genes=which(table_EB_zigFit_contrasts$adj.P.Val<0.05)
with(table_EB_zigFit_contrasts, symbols(logFC, -log10(adj.P.Val), pch=20, 
                            circles = radius_phylum_node,inches=0.1, fg="black",
                            bg=ifelse(table_EB_zigFit_contrasts$adj.P.Val<.06,"red","gray"),
                            main="Volcano plot"))

#
## Here's another example of a volcano plot
volcanoplot(EB_zigFit_contrasts, highlight = 8, names = row.names(table_EB_zigFit_contrasts))

