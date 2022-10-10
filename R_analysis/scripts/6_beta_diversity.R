# Load previous work
source("Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/scripts/5_intro_to_ggplot.R")

#### RECAP OF PHYLOSEQ OBJECTS ####
# Let's revist what we have done so far, and what we are going to use for β-diversity.
# Here are the phyloseq objects we have:
microbiome.ps # not normalized; 16S ASV counts, tree, rep-seqs, taxonomy, tree, metadata
amr.ps # not normalized; ARG counts, taxonomy, metadata

microbiome_css.ps # CSS normalized; 16S ASV counts, tree, rep-seqs, taxonomy, tree, metadata
amr_css.ps # CCS normalized; ARG counts, taxonomy, metadata

microbiome_ccs_phylum.ps # CSS normalized, agglomerated at Phylum, not normalized within sample (RA)
amr_ccs_phylum.ps # CSS normalized, agglomerated at class, not normalized within sample (RA)

# Which should we use moving forward to analyze beta-diversity?
# **We could use the agglomerated counts to make RA plots for Phylum and Class
# But for ordination we need the non-agglomerated counts

#### ORDINATION ####

# Step 1 for both of these is the calculate a distance matrix
# Lots of options, but we are going to use Bray-Curtis and the vegan package's 'vegdist' function
microbiome_bray.dist <- vegdist(t(otu_table(microbiome_css.ps)), method = "bray")
amr_bray.dist <- vegdist(t(otu_table(amr_css.ps)), method = "bray")

# take a peek
microbiome_bray.dist
amr_bray.dist

# Step 2 is to ordinate the distance matrices using phyloseq's 'ordinate'
microbiome_bray.ord <- ordinate(microbiome_css.ps, method = "NMDS", distance = microbiome_bray.dist)
amr_bray.ord <- ordinate(amr_css.ps, method = "NMDS", distance = amr_bray.dist)
# Note: we can ordinate with a bunch of different methods:
# PCoA, CCA, RDA, MDS, etc.

# Step 3 is to plot the NMDS using phyloseq's 'plot_ordination'

#### 16S
plot_ordination(microbiome_css.ps, microbiome_bray.ord, type = "samples", color = "Sample_type")
# WWTP communities are SO different its tough to visualize
# We can use facet_wrap to help
# but still saving the plot for stats section
microbiome_NMDS_plot <- plot_ordination(microbiome_css.ps, microbiome_bray.ord, type = "samples", color = "Sample_type") +
  theme_bw() +
  stat_ellipse(lty = 2) +
  geom_point(size = 4, shape = 18) +
  scale_color_manual(values = c("dodgerblue3","orange3","mediumorchid3","seagreen3"))


# IMPORTANT POINT FOR INTEPRETATION VS. VISUALIZATION
plot_ordination(microbiome_css.ps, microbiome_bray.ord, type = "samples", color = "Sample_type") +
  facet_wrap(~Sample_type)

plot_ordination(microbiome_css.ps, microbiome_bray.ord, type = "samples", color = "Sample_type") +
  facet_wrap(~Sample_type, scales = "free")

# AMR-TE
plot_ordination(amr_css.ps, amr_bray.ord, type = "samples", color = "Sample_type")
# much more normal looking; distance between WWTP and others not as large (still large, though!)

# Moving forward, we're just going to look at the AMR-TE for NMDS

# Remember, we have the many plotting options for this "phyloseq plot"
# For example, we can add a 95% confidence ellipses around our sample groups,
# but this doesn't work well if there are no clusters
# We can also make the points larger and change the shape and color
amr_NMDS_plot <- plot_ordination(amr_css.ps, amr_bray.ord, type = "samples", color = "Sample_type") +
  theme_bw() +
  stat_ellipse(lty = 2) +
  geom_point(size = 4, shape = 18) +
  scale_color_manual(values = c("dodgerblue3","orange3","mediumorchid3","seagreen3"))
amr_NMDS_plot

# another example with facetting
plot_ordination(amr_css.ps, amr_bray.ord, type = "samples", color = "Sample_type") +
  theme_bw() + facet_wrap(~Sample_type) +
  stat_ellipse(lty = 2) +
  geom_point(size = 4, shape = 18) +
  scale_color_manual(values = c("dodgerblue3","orange3","mediumorchid3","seagreen3"))


#### HIERARCHAL CLUSTERING ####
# For this we start with the distance matrix we created above
# Step 1 is to perform hierarchal clustering, which can be done with different methods
# For this we're going to use Ward's agglomeration method (see the help for 'hclust' for other options)
microbiome.hclust <- hclust(microbiome_bray.dist, method = "ward.D2")
amr.hclust <- hclust(amr_bray.dist, method = "ward.D2")
plot(microbiome.hclust)
plot(amr.hclust)

# plot is fine, but if we want to make it nicer using ggplot
ggplot(amr.hclust, aes(x=x,y=y)) # not possible to use ggplot to plot hclust objects

# So, Step2 is to convert to a dendrogram and use the packaged ggdendro to plot including metadata info
microbiome.dendro <- as.dendrogram(microbiome.hclust)
amr.dendro <- as.dendrogram(amr.hclust)
plot(microbiome.dendro) # no difference, just lost the axes etc.

# Step 3 is to get the plotting data (co-ordinates for the lines etc.) and 
# the labels (i.e., metadata, which at this point is just sample name) from our dendro object
# so we can add additional metadata details to the labels for meaningful plots
# grabbing the data
microbiome.dendro.data <- dendro_data(microbiome.dendro, type = "rectangle")
amr.dendro.data <- dendro_data(amr.dendro, type = "rectangle")
# depending on what kind of dendro you want to plot you can choose 'rectangle' or 'triangle'

# Step 4 is amending our metadata to the dendrogram labels

#16S
# let's look at the current labels(we've got SampleID there and it's named 'label')
microbiome.dendro.data$labels
# convert metadata to tibble data frame
microbiome_metadata_for_dendro <- as_tibble(microbiome_css.ps@sam_data)
# use left_join to add metadata to the dendro.data$labels
microbiome.dendro.data$labels <- microbiome.dendro.data$labels %>%
  left_join(microbiome_metadata_for_dendro, by = c("label" = "SampleID"))
# look at the new labels
microbiome.dendro.data$labels # got all the metadata attached!

# AMR-TE
# convert metadata to tibble data frame
amr_metadata_for_dendro <- as_tibble(sample_data(sample_data(amr_css.ps)))
# use left_join to add metadata to the dendro.data$labels
amr.dendro.data$labels <- amr.dendro.data$labels %>%
  left_join(amr_metadata_for_dendro, by = c("label" = "SampleID"))
# look at the new labels
amr.dendro.data$labels # got all the metadata attached!

# PLOT THE DENDROGRAMS
# now we can use ggplot and can add colours, shapes, etc based on our metadata
microbiome_dendro_plot <- ggplot(microbiome.dendro.data$segments) +
  theme_minimal() +
  labs(title = "Microbiome Hierarchal Clustering", y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = microbiome.dendro.data$labels, aes(x=x,y=y, colour = Sample_type), size = 12, shape = 15) +
  theme(plot.title = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
microbiome_dendro_plot

# now we can use ggplot and can add colours, shapes, etc based on our metadata
amr_dendro_plot <- ggplot(amr.dendro.data$segments) +
  theme_minimal() +
  labs(title = "Resistome Hierarchal Clustering", y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = microbiome.dendro.data$labels, aes(x=x,y=y, colour = Sample_type), size = 12, shape = 15) +
  theme(plot.title = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
amr_dendro_plot

#### RELATIVE ABUNDANCE & RELATED FIGURES ####

# One of the easiest way we can begin to compare the microbiome and resistome composition is by plotting relative abundance
# Convert OTU abundances to relative abundances (I prefer making it a % opposed to proportion, but up to you)
microbiome_ra.ps <- transform_sample_counts(microbiome_css.ps, function(x) {x/sum(x)}*100)
amr_ra.ps <- transform_sample_counts(amr_css.ps, function(x) {x/sum(x)}*100)

#if you wanted to double check, all sample sums (total ASVs within one sample) should equal 100
plot(sample_sums(microbiome_ra.ps), type = "h", ylab = "total ASVs", xlab= "samples")
plot(sample_sums(microbiome_ra.ps), type = "h", ylab = "total ASVs", xlab= "samples")

# agglomerate at the phylum and class level
microbiome_phylum_ra.ps <- tax_glom(microbiome_ra.ps, taxrank = "Phylum", NArm = F)
microbiome_phylum_ra.ps # 37 phyla (same as before)
amr_class_ra.ps <- tax_glom(amr_ra.ps, taxrank = "class")
amr_class_ra.ps # 48 classes (same as before)

# plot this (notice the y-axis as compared to last time we plotted these)
plot_bar(microbiome_phylum_ra.ps, fill="Phylum")
plot_bar(amr_class_ra.ps, fill="class") # this is tough to see with the legend this big

# you could use these figures, but in order to make them nice using ggplot
# we need to melt the data into "long" format data frame
microbiome_phylum_ra.melt <- psmelt(microbiome_phylum_ra.ps)
amr_class_ra.melt <- psmelt(amr_class_ra.ps)

# we can see we still have 37 phylum and
unique(microbiome_phylum_ra.melt$Phylum) # there are our 37 phyla
unique(amr_class_ra.melt$class) # and our 48 classes

# now we can use ggplot to plot these
# can plot the mean RA by sample type
ggplot(microbiome_phylum_ra.melt, aes(x= Sample_type, y= Abundance, fill = Phylum)) +
  theme_bw() +
  geom_bar(stat = "summary", colour = "black")

ggplot(amr_class_ra.melt, aes(x= Sample_type, y= Abundance, fill = class)) +
  theme_bw() +
  geom_bar(stat = "summary", colour = "black")
# hard to see the plot with this many taxa
# So lets change the names of phyla to "low abundance phyla" for features present at a relative
# proportion less than 0.1 across all samples
microbiome_phylum_ra.melt <- microbiome_phylum_ra.melt %>%
  group_by(Phylum) %>%
  mutate(mean_phylum_rel_abundance = mean(Abundance))

microbiome_phylum_ra.melt$Phylum <- as.character(microbiome_phylum_ra.melt$Phylum)
microbiome_phylum_ra.melt$mean_phylum_rel_abundance<- as.numeric(microbiome_phylum_ra.melt$mean_phylum_rel_abundance)

microbiome_phylum_ra.melt$Phylum[microbiome_phylum_ra.melt$mean_phylum_rel_abundance < 0.1] <- "Low abundance phyla"
unique(microbiome_phylum_ra.melt$Phylum) # down to 21 phyla from 37

# Same thing now for classes of AMR
amr_class_ra.melt <- amr_class_ra.melt %>%
  group_by(class) %>%
  mutate(mean_class_rel_abundance = mean(Abundance))

amr_class_ra.melt$class<- as.character(amr_class_ra.melt$class)
amr_class_ra.melt$mean_class_rel_abundance<- as.numeric(amr_class_ra.melt$mean_class_rel_abundance)

amr_class_ra.melt$class[amr_class_ra.melt$mean_class_rel_abundance < 0.1] <- "Low abundance classes"
unique(amr_class_ra.melt$class) # down to 27 classes from 48

# re-plot them with low abundance phyla combined
microbiome_phylaRA_sampleType_plot <- ggplot(microbiome_phylum_ra.melt, aes(x= Sample_type, y= Abundance, fill = Phylum)) +
  theme_bw() +
  geom_bar(stat = "summary", colour = "black")
microbiome_phylaRA_sampleType_plot

amr_classRA_sampleType_plot <- ggplot(amr_class_ra.melt, aes(x= Sample_type, y= Abundance, fill = class)) +
  theme_bw() +
  geom_bar(stat = "summary", colour = "black")
amr_classRA_sampleType_plot

#### COMBINING RA WITH OUR DENDROGRAMS ####
# For this, we want barplots for each invidivdual sample, and in the order of they are
# within the dendrogram

# use the dendrogram label data to position the samples in the right order
# make object of microbiome sample names
microbiome_sample_names <- microbiome.dendro.data$labels$label
# plot the RA barplot by individual samples
microbiome_ra_dendro_plot <- ggplot(microbiome_phylum_ra.melt, aes(x= SampleID, y= Abundance, fill = Phylum)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = microbiome_sample_names) +
  theme(plot.title = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
microbiome_ra_dendro_plot
# combine this plot with the dendrogram we made earlier
plot_grid(microbiome_dendro_plot, microbiome_ra_dendro_plot, align = "v", ncol = 1)

# make object of amr sample names
amr_sample_names <- amr.dendro.data$labels$label
# plot the RA barplot
amr_ra_dendro_plot = ggplot(amr_class_ra.melt, aes(x= SampleID, y= Abundance, fill = class)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = amr_sample_names) +
  theme(plot.title = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
amr_ra_dendro_plot

# combine this plot with the dendrogram we made earlier
plot_grid(amr_dendro_plot, amr_ra_dendro_plot, align = "v", ncol = 1)


#### TARGETTING SPECIFIC TAXA OF INTEREST ####

# Sometimes we might be interested in specific taxa
# For example, in our 16S data we can see that the RA of phylum Bacteroidota is different between samples
# or you may be interested in a particular taxa for your research question
# let's pretend we're interested in Bacteroidota for 16S and Tetracyclines for AMR
microbiome_phylaRA_sampleType_plot
amr_classRA_sampleType_plot

#Bacteroidota
bacteroidota_phylum <- microbiome_phylum_ra.melt[which(microbiome_phylum_ra.melt$Phylum=="Bacteroidota"),]
bacteroidota_phylum$Phylum # all we've got is Bacteroidota

bacteroidota_sampletype_plot <- ggplot(bacteroidota_phylum, aes(x= Sample_type, y= Abundance, fill = Sample_type)) +
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "Bacteroidota") +
  geom_boxplot() +
  geom_point()
bacteroidota_sampletype_plot

# Tetracyclines
tetracycline_class <- amr_class_ra.melt[which(amr_class_ra.melt$class=="Tetracyclines"),]
tetracycline_class$class # good, we've only got that class of ARGs

tetracycline_sampletype_plot <- ggplot(tetracycline_class, aes(x= Sample_type, y= Abundance, fill = Sample_type)) +
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "Genes confering resistance to Tetracyclines") +
  geom_boxplot() +
  geom_point()
tetracycline_sampletype_plot
