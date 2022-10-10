#load previous work
source("Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/scripts/2_summary_statistics.R")

#### CALCULATING ALPHA-DIVERSITY METRICS ####

# Alpha-diversity indices are a common measurement to summarize the composition of the microbiome and resistome.
# "Richness" or "Observed" simply describes the number of unique taxa identified in a sample
# On the other hand, we use "Shannon's index" to describe both richness and "evenness", or
# distribution of counts among taxa in a sample

# Find more infomation about diversity indices here: http://www.jmb.or.kr/submission/Journal/027/JMB027-12-02_FDOC_2.pdf

# calculating richness and Shannon
microbiome_alpha <- estimate_richness(microbiome.ps, measures = c("Observed","Shannon"))
microbiome_alpha

amr_alpha <- estimate_richness(amr.ps, measures = c("Observed","Shannon"))
amr_alpha

# to make meaningful comparisons we need to include the metadata
# convert metadata into a data.frame
microbiome_metadata <- as(sample_data(microbiome.ps),"data.frame")
amr_metadata <- as(sample_data(amr.ps),"data.frame")
# combine the metadata and alpha diversity values
microbiome_alpha_meta <- cbind(microbiome_metadata, microbiome_alpha)
microbiome_alpha_meta
amr_alpha_meta <- cbind(amr_metadata, amr_alpha)
amr_alpha_meta

#### PLOTTING RICHNESS ####
microbiome_richness_boxplot <- ggplot(microbiome_alpha_meta, aes(x= Sample_type, y= Observed, fill = Sample_type)) +
  theme_bw() + labs(y= "Observed ASVs") +
  geom_boxplot() +
  geom_point()
microbiome_richness_boxplot

amr_richness_boxplot <- ggplot(amr_alpha_meta, aes(x= Sample_type, y= Observed, fill = Sample_type)) +
  theme_bw() + labs(y= "Observed ARGs") +
  geom_boxplot() +
  geom_point()
amr_richness_boxplot

#### PLOTTING DIVERSITY ####
microbiome_diversity_boxplot <- ggplot(microbiome_alpha_meta, aes(x= Sample_type, y= Shannon, fill = Sample_type)) +
  theme_bw() + labs(y= "Observed ASVs") +
  geom_boxplot() +
  geom_point()
microbiome_diversity_boxplot

amr_diversity_boxplot <- ggplot(amr_alpha_meta, aes(x= Sample_type, y= Shannon, fill = Sample_type)) +
  theme_bw() + labs(y= "Observed ARGs") +
  geom_boxplot() +
  geom_point()
amr_diversity_boxplot
