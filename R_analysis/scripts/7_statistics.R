####load previous scripts ####
source("Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/scripts/6_beta_diversity.R")


#### STATISTICS ####

#### COMPARING MEANS - NON-PARAMETRIC ####
# Microbiome/resistome is not normally distributed
# We need to use non-parametric tests that do not assume normal distribution

### pairwise wilcoxon rank sum is used for anything comparing means of variable

### let's repeat the test on sequencing depth between sample types
# 16S
sequencing_depth_16S_boxplot
pairwise.wilcox.test(sample_data(microbiome.ps)$Raw_paired_reads, sample_data(microbiome.ps)$Sample_type)
# No signficant differences using the holm adjustment for multiple comparisons, which we can switch
pairwise.wilcox.test(sample_data(microbiome.ps)$Raw_paired_reads, sample_data(microbiome.ps)$Sample_type, p.adjust.method = "BH")
# p-values different, but no significant differences

# AMR-TE
sequencing_depth_amr_boxplot
pairwise.wilcox.test(sample_data(amr.ps)$Raw_paired_reads, sample_data(amr.ps)$Sample_type, p.adjust.method = "BH")
# No signficant differences

#### STATS ON ALPHA DIVERSITY ####
# 16S
# We can also use the pairwise wilcoxon for testing alpha diversity metrics
microbiome_richness_boxplot
pairwise.wilcox.test(microbiome_alpha_meta$Observed, microbiome_alpha_meta$Sample_type, p.adjust.method = "BH")
# despite looking different, no signficant differencces in microbial richness
microbiome_diversity_boxplot
pairwise.wilcox.test(microbiome_alpha_meta$Shannon, microbiome_alpha_meta$Sample_type, p.adjust.method = "BH")
# again, no significant differences

#AMR-TE
# We can also use the pairwise wilcoxon for testing alpha diversity metrics
amr_richness_boxplot
pairwise.wilcox.test(amr_alpha_meta$Observed, amr_alpha_meta$Sample_type, p.adjust.method = "BH")
# all sig. diff. from all
amr_diversity_boxplot
pairwise.wilcox.test(amr_alpha_meta$Shannon, amr_alpha_meta$Sample_type, p.adjust.method = "BH")
# all significantly different from each other, except swine and WWTP not different

# We can also use this to test for differences in individual taxa
# Remember we looked at Bacteroidota and Tetracyclines
bacter

#### PERMUTATIONAL MANOVA (adonis) ####
# We use PERMANOVA for beta-diversity
# 16S
microbiome_NMDS_plot
microbiome

# We'll need to convert our data to a data.frame object
# There is probably a better way to extract the data than using
# all of the nested functions you'll see below, but this works for us!
microbiome_df = as.data.frame(as(sample_data(microbiome_css.ps),"matrix"))

# To run the PERMANOVA we need the dist matrix we calculated previously
microbiome_bray.dist
# create the PERMANOVA model
adonis_microbiome_sampleType <- pairwise.adonis2(microbiome_bray.dist ~ Sample_type, microbiome_df)
# view the results
adonis_microbiome_sampleType # significant, But...

# Check assumption of similar multivariate spread among the treatment groups.
# PERMANOVA works with the assumption that dispersion of the data in your 
# samples is the ~same among each other, so before running PERMANOVA,
# you must run the betadisper->permutest to know if the dispersions are the 
# same. For that to be true, the permutest has to have a 
# non-significant p-value. Knowing the previous, 
# then you can run the PERMANOVA test, otherwise your interpretations will be wrong”
permutest(betadisper(microbiome_bray.dist, microbiome_df$Sample_type), pairwise = TRUE)
# we have some signficant differences, meaning we cannot rule out the significance of
# differences in dispersion, however we can still draw conclusions

# same thing with the AMR-TE
amr_NMDS_plot
amr_df <- as.data.frame(as(sample_data(amr_css.ps),"matrix"))
adonis_amr_sampleType <- pairwise.adonis2(amr_bray.dist ~ Sample_type, amr_df)
adonis_amr_sampleType # all from all
permutest(betadisper(amr_bray.dist, amr_df$Sample_type), pairwise = TRUE)
# no significant differences in dispersions of variance!

