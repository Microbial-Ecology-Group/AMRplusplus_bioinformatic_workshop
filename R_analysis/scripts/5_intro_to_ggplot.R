# Load previous work
source("Desktop/AMRplusplus_bioinformatic_workshop/R_analysis/scripts/4_count_normalization.R")

# recall the boxplot for amr sequencing depth?
sequencing_depth_amr_boxplot


#### GGPLOT LAYERS ####
# We can create the initial plot by specifying the data that we'll be using and the x and y variables (columns)
ggplot(sample_data(amr.ps), aes(x = Sample_type , y = Raw_paired_reads))

# A blank ggplot is drawn. Even though the x and y are specified, there are no points or lines in it.
# This is because, ggplot doesn’t assume that you meant a scatterplot or a line chart to be drawn.
# We've only told ggplot what dataset to use and what columns should be used for X and Y axis, but we 
# haven’t explicitly asked it to draw any points.

# Also note that aes() function is used to specify the X and Y axes. 
# That’s because, any information that is part of the source dataframe has to be specified inside the aes() function.

# To create a boxplot, we can add a "geom" layer, geom_boxplot()
# Other options include (and many others)
# geom_point()
# geom_line()
# geom_bar()

# So, let's reconstruct our sequencing depth boxplot
# Here's an example of a bare bones boxplot
ggplot(sample_data(amr.ps), aes(x = Sample_type , y = Raw_paired_reads)) + 
  geom_boxplot()

# Or we can turn this figure into a graph with points
ggplot(sample_data(amr.ps), aes(x = Sample_type , y = Raw_paired_reads)) + 
  geom_point()

# Or we can do both
ggplot(sample_data(amr.ps), aes(x = Sample_type , y = Raw_paired_reads)) + 
  geom_boxplot() +
  geom_point()

# We can add color based on another variable
ggplot(sample_data(amr.ps), aes(x = Sample_type , y = Raw_paired_reads, color = Sample_type)) + 
  geom_boxplot() +
  geom_point()

# We can add fill color based on another variable
ggplot(sample_data(amr.ps), aes(x = Sample_type , y = Raw_paired_reads, fill = Sample_type)) + 
  geom_boxplot() +
  geom_point()

# Specifying the aesthetics "aes()" in the plot vs. in the layers
# Aesthetic mappings can be supplied in the initial ggplot() call, in individual layers, or in some combination of both. 
# All of these calls create the same plot specification:
ggplot(sample_data(amr.ps), aes(Sample_type, Raw_paired_reads, fill = Sample_type)) + 
  geom_boxplot()
ggplot(sample_data(amr.ps), aes(Sample_type, Raw_paired_reads)) + 
  geom_boxplot(aes(fill = Sample_type))
ggplot(sample_data(amr.ps), aes(Sample_type)) + 
  geom_boxplot(aes(y = Raw_paired_reads, fill = Sample_type))
ggplot(sample_data(amr.ps)) + 
  geom_boxplot(aes(Sample_type, Raw_paired_reads, fill = Sample_type))

# Each of the labels has different arguments you can use to change the aesthetics.
ggplot(sample_data(amr.ps), aes(Sample_type, Raw_paired_reads, fill = Sample_type)) + 
  geom_boxplot(colour= "green")

# Or we can color based on unique factors
ggplot(sample_data(amr.ps), aes(Sample_type, Raw_paired_reads, fill = Sample_type)) + 
  geom_boxplot(aes(colour= Sample_type))

# There are many options for layers to add, so be sure to search online for your specific needs
ggplot(sample_data(amr.ps), aes(Sample_type, Raw_paired_reads, fill = Sample_type)) + 
  geom_boxplot(colour = "black") +
  geom_jitter(width = 0.1)

#### GGPLOT LAYERS ####
# We can also add layers with labels for the graph using ggtitle(), xlab(), ylab(), and labs
ggplot(sample_data(amr.ps), aes(x = Sample_type , y = Raw_paired_reads, fill = Sample_type)) + 
  geom_boxplot() +
  geom_point() +
  ggtitle("Metagenomic sequencing results")

# We can also change the title,  x and y axis labels with labs()
ggplot(sample_data(amr.ps), aes(x = Sample_type , y = Raw_paired_reads, fill = Sample_type)) + 
  geom_boxplot() +
  geom_point() +
  labs(x = "Sample type", y = "raw paired reads", title = "Metagenomic sequencing results")

# Or we can use xlab() and ylab()
ggplot(sample_data(amr.ps), aes(x = Sample_type , y = Raw_paired_reads, fill = Sample_type)) + 
  geom_boxplot() +
  geom_point() +
  ggtitle("Metagenomic sequencing results") +
  ylab("raw paired reads") +
  xlab("Sample type")

# You can add a caption
ggplot(sample_data(amr.ps), aes(x = Sample_type , y = Raw_paired_reads, fill = Sample_type)) + 
  geom_boxplot() +
  geom_point() +
  labs(caption = "(based on data from...)")

# Or a tag, for panel figures
ggplot(sample_data(amr.ps), aes(x = Sample_type , y = Raw_paired_reads, fill = Sample_type)) + 
  geom_boxplot() +
  geom_point() +
  labs(tag = "A")

#### GGPLOT THEMES ####
# ggplot has a lot of ways to alter the color schemes and the overall aesthetic of the graph
# we like using "theme()" as an easy way to try other color formatting

# You can alter each specific part of the theme, or use predefined themes such as "theme_classic()" shown below
# Again, there are so many options for what you can modify in the theme(), so look up what you need online and you'll probably find it

ggplot(sample_data(amr.ps), aes(x= Sample_type, y= Raw_paired_reads, fill = Sample_type)) +
  labs(title = "AMR-TE sequencing depth", y= "No. raw paired reads", x= "Sample type") +
  geom_boxplot() +
  geom_point() +
  theme(plot.title = element_text(size = 24),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# theme_classic
ggplot(sample_data(amr.ps), aes(x= Sample_type, y= Raw_paired_reads, fill = Sample_type)) +
  theme_classic() +
  labs(title = "AMR-TE sequencing depth", y= "No. raw paired reads", x= "Sample type") +
  geom_boxplot() +
  geom_point()

# theme_bw
ggplot(sample_data(amr.ps), aes(x= Sample_type, y= Raw_paired_reads, fill = Sample_type)) +
  theme_bw() +
  labs(title = "AMR-TE sequencing depth", y= "No. raw paired reads", x= "Sample type") +
  geom_boxplot() +
  geom_point()

# theme_minimal
ggplot(sample_data(amr.ps), aes(x= Sample_type, y= Raw_paired_reads, fill = Sample_type)) +
  theme_minimal() +
  labs(title = "AMR-TE sequencing depth", y= "No. raw paired reads", x= "Sample type") +
  geom_boxplot() +
  geom_point()

# theme_dark
ggplot(sample_data(amr.ps), aes(x= Sample_type, y= Raw_paired_reads, fill = Sample_type)) +
  theme_dark() +
  labs(title = "AMR-TE sequencing depth", y= "No. raw paired reads", x= "Sample type") +
  geom_boxplot() +
  geom_point()

# You can also override theme elements by combining the two options
# regular theme_bw()
ggplot(sample_data(amr.ps), aes(x= Sample_type, y= Raw_paired_reads, fill = Sample_type)) +
  theme_bw() +
  labs(title = "AMR-TE sequencing depth", y= "No. raw paired reads", x= "Sample type") +
  geom_boxplot() +
  geom_point()

# customized theme_bw()
ggplot(sample_data(amr.ps), aes(x= Sample_type, y= Raw_paired_reads, fill = Sample_type)) +
  theme_bw() +
  labs(title = "AMR-TE sequencing depth", y= "No. raw paired reads", x= "Sample type") +
  geom_boxplot() +
  geom_point() +
  theme(plot.title = element_text(size = 24),
        panel.border = element_rect(colour = "black", size = 0.75),
        panel.grid.major.x = element_blank(),
        legend.key.size = unit(2,"lines"),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#### GGPLOT FACETS ####
# we can split the data into different figures using facets
ggplot(sample_data(amr.ps), aes(x= ID, y= Raw_paired_reads, fill = Sample_type)) +
  theme_bw() + facet_wrap(~Sample_type) +
  labs(title = "AMR-TE sequencing depth", y= "No. raw paired reads", x= "Sample type") +
  geom_bar(stat = "summary") +
  geom_point() +
  theme(plot.title = element_text(size = 24),
        panel.border = element_rect(colour = "black", size = 0.75),
        panel.grid.major.x = element_blank(),
        legend.key.size = unit(2,"lines"),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# we can also make the scales free free to change depending on the facet
ggplot(sample_data(amr.ps), aes(x= ID, y= Raw_paired_reads, fill = Sample_type)) +
  theme_bw() + facet_wrap(~Sample_type, scales = "free_x") +
  labs(title = "AMR-TE sequencing depth", y= "No. raw paired reads", x= "Sample type") +
  geom_bar(stat = "summary") +
  geom_point() +
  theme(plot.title = element_text(size = 24),
        panel.border = element_rect(colour = "black", size = 0.75),
        panel.grid.major.x = element_blank(),
        legend.key.size = unit(2,"lines"),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
  