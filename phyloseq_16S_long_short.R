#Load libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(readxl)

# Set the path to the BIOM file
biom_file <- "/Microbiome Stage 1/Prelimstudy_write_up/Export_merged/table-with-taxonomy.biom"

# Import the BIOM file into a phyloseq object
psPrelim <- import_biom(biom_file)
psPrelim

#Check tax labels in BIOM
colnames(tax_table(psPrelim))

#Check sample labels in BIOM
colnames(otu_table(psPrelim))

#rename tax labels
colnames(tax_table(psPrelim)) <- c("Kingdom", "Phylum", "Class", 
                             "Order", "Family", "Genus", "Species")

# Import sample metadata
mapPrelim <- read_xlsx("Desktop/Microbiome Stage 1/Prelimstudy_write_up/Manifest_All_prelim_study_forfeturetable.xlsx")
mapPrelim


#convert year to a factor
mapPrelim[, 'Year'] <- lapply(mapPrelim[, 'Year'], factor)
mapPrelim

mapPrelim <- sample_data(mapPrelim)
mapPrelim


#Assign row names to be Sample ID's
row.names(mapPrelim) <- mapPrelim$id
mapPrelim

sample_names(mapPrelim)

sample_names(psPrelim)


# Merge data object with sample metadata
data_mergePrelim <- merge_phyloseq(psPrelim, mapPrelim)
data_mergePrelim

sample_names(data_mergePrelim)

colnames(tax_table(data_mergePrelim))


OTU_table_Prelim <- otu_table(data_mergePrelim)
head(OTU_table_Prelim)

#Write csv
write.csv(OTU_table_Prelim, file = "Phyloseq_OTU_tablePrelim.csv")


#First look at the distribution of read counts from our samples
# Make a data frame with a column for the read counts of each sample

sample_sum_df <- data.frame(sum = sample_sums(data_mergePrelim))
sample_sum_df

# Save the data frame as a CSV file
write.csv(sample_sum_df, "/Writing/Prelim study/Extra figures/Seq_deph_all.csv", row.names = TRUE)

number_of_asvs <- ntaxa(physeq)

# Print the number of ASVs
print(number_of_asvs)

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  theme_bw() +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())


read_sums_Merged_data_from_phyloseq <- read_excel("/Prelimstudy_write_up/Illumina_PacBio_merge_for_anal/read_sums_Merged_data_from_phyloseq.xlsx")
read_sums_Merged_data_from_phyloseq

# Histogram for PacBio
His1 <- ggplot(data = subset(read_sums_Merged_data_from_phyloseq, Platform == "PacBio"), aes(x = sum)) +
  geom_histogram(binwidth = 500, fill = "blue", color = "black") +
  theme_bw() +
  labs(title = "Read counts PacBio", x = "Sum", y = "Frequency") +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

His1

# Histogram for Other Platform
His2 <- ggplot(data = subset(read_sums_Merged_data_from_phyloseq, Platform != "PacBio"), aes(x = sum)) +
  geom_histogram(binwidth = 500, fill = "purple", color = "black") +
  theme_bw() +
  labs(title = "Read counts Illumina", x = "Sum", y = "Frequency") +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

His2

library(gridExtra)

# Arrange the histograms together
grid.arrange(His1, His2, ncol = 2)

# mean, max and min of sample read counts
smin <- min(sample_sums(data_mergePrelim))

smean <- mean(sample_sums(data_mergePrelim))
smax <- max(sample_sums(data_mergePrelim))

smean
smax
# melt to long format (for ggploting) 
# prune out phyla below 1% in each sample


BIOMdata_merge_Phylum <- data_mergePrelim %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.001) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# Set colors for plotting
phylum_colors <- c(
  "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c","#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",
  "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#5F7FC7", "orange","#DA5724", "#508578",
  "#CD9BCD","#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#CBD588","darkgray", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "blue",
  "orange","green","cornsilk","darkblue","darkturquoise","brown4","deeppink4", "aliceblue", "cadetblue","darkgreen"
)

BIOMdata_merge_Phylum

# Plot 
ggplot(BIOMdata_merge_Phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(
    breaks = unique(BIOMdata_merge_Phylum$Platform),  # Set breaks to unique platform values
    labels = unique(BIOMdata_merge_Phylum$Platform),  # Set labels to unique platform values
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 0.001%) \n") +
  ggtitle("Phylum Composition of Soil Communities by Sample Number") 


ggplot(BIOMdata_merge_Phylum, aes(x = factor(Platform), y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(
    breaks = unique(BIOMdata_merge_Phylum$Platform),  # Set breaks to unique platform values
    labels = unique(BIOMdata_merge_Phylum$Platform),  # Set labels to unique platform values
    drop = FALSE
  ) +
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 0.001%) \n") +
  ggtitle("Phylum Composition of Soil Communities by Platform")


ggplot(BIOMdata_merge_Phylum, aes(x = factor(Platform), y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(
    breaks = unique(BIOMdata_merge_Phylum$Platform),
    labels = unique(BIOMdata_merge_Phylum$Platform),
    drop = FALSE
  ) +
  theme_bw() +
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, order = -1, ncol = 2)) +  # Specify 2 columns
  ylab("Relative Abundance (Phyla > 0.001%) \n") +
  ggtitle("Phylum Composition of Soil Communities by Platform") +
  scale_y_continuous(breaks = seq(0, 110, by = 10)) + # Add more numbers to y-axis
theme(axis.text.x = element_text(size = 16, face='bold'))+
  theme(axis.title=element_text(size=16, face='bold'),
        axis.text=element_text(size=16, face='bold'),
        strip.text=element_text(size=16, face='bold'),
        strip.background=element_rect(fill=NA, color='black', size=0.5),
        panel.border=element_rect(fill=NA, color='black')) + 
  theme(legend.position = "right", 
        plot.title= element_text(size = 16, hjust = 0.5, face='bold')) 


##########Other_PLOTS###########
# Set colors for plotting
phylum_colors <- c(
  "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c","#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",
  "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#5F7FC7", "orange","#DA5724", "#508578",
  "#CD9BCD","#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#CBD588","darkgray", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "blue",
  "orange","green","cornsilk","darkblue","darkturquoise","brown4","deeppink4", "aliceblue", "cadetblue","darkgreen"
)

#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(data_mergePrelim))
standf = function(x, t=total) round(t * (x / sum(x)))
data_mergePrelim = transform_sample_counts(data_mergePrelim_filtered, standf)
data_mergePrelim

BIOMdata_merge_Phylum2 <- data_mergePrelim %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.001) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum


ggplot(BIOMdata_merge_Phylum2, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(
    breaks = unique(BIOMdata_merge_Phylum$Platform),  # Set breaks to unique platform values
    labels = unique(BIOMdata_merge_Phylum$Platform),  # Set labels to unique platform values
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 0.001%) \n") +
  ggtitle("Phylum Composition of Soil Communities by Sample Number") 

sample_names(data_mergePrelim)
sample_variables(data_mergePrelim)

Prelim_Seq_platform <- merge_samples(data_mergePrelim, "Platform")


BIOMdata_merge_Phylum2 <- Prelim_Seq_platform %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.001) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum



ggplot(BIOMdata_merge_Phylum2, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(
    breaks = unique(BIOMdata_merge_Phylum$Platform),  # Set breaks to unique platform values
    labels = unique(BIOMdata_merge_Phylum$Platform),  # Set labels to unique platform values
    drop = FALSE
  ) +
  theme_bw() +
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, order = -1, ncol = 2)) +  # Specify 2 columns
  ylab("Relative Abundance (Phyla > 0.001%) \n") +
  ggtitle("Phylum Composition of Soil Communities by Sequencing Platform") +
  theme(axis.text.x = element_text(size = 14, face='bold'))+
  theme(axis.title=element_text(size=14, face='bold'),
        axis.text=element_text(size=14, face='bold'),
        strip.text=element_text(size=14, face='bold'),
        strip.background=element_rect(fill=NA, color='black', size=0.5),
        panel.border=element_rect(fill=NA, color='black')) + 
  theme(legend.position = "right", 
        plot.title= element_text(size = 14, hjust = 0.5, face='bold')) 
  


#######################
# Ordinate
BIOMdatamerge_pcoa <- ordinate(
  physeq = data_mergePrelim, 
  method = "PCoA", 
  distance = "bray"
)

#PCoA by year
# Plot 
p1a2 = plot_ordination(
  physeq = data_mergePrelim,
  ordination = BIOMdatamerge_pcoa,
  color = "Year",
  type="samples",
  title = "PCoA of bacterial Communities"
) +
  theme_bw() +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                "#CAB2D6", "#6A3D9A", "white", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD")) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Platform), type = "norm", level = 0.95, color = "#999999", size = 0.5) +
  labs(x = "PCoA axis 1 (15.4%)",
       y = "PCoA axis 1 (5.7%)") + 
  theme(axis.text.x = element_text(size = 16, face='bold'))+
  theme(axis.title=element_text(size=16, face='bold'),
        axis.text=element_text(size=16, face='bold'),
        strip.text=element_text(size=16, face='bold')) + 
  theme(legend.position = "right", 
        plot.title= element_text(size = 16, hjust = 0.5, face='bold')) +
  

p1a2




Tree <- "/Microbiome Stage 1/Prelimstudy_write_up/Export_merged/tree.nwk" 

# Read the tree file
library("ape")
tree <- read.tree(Tree)

Tree2 <- "/Microbiome Stage 1/Prelimstudy_write_up/Export_merged/tree_rooted.nwk" 

# Read the tree file
library("ape")
tree2 <- read.tree(Tree2)
########Rebuild phyloseq data from scratch using all the simulated data components we just generated###########
  
physeq2 = phyloseq(OTU_table_Prelim, psPrelim, mapPrelim, tree)
physeq2

physeq3 = phyloseq(OTU_table_Prelim, psPrelim, mapPrelim, tree2)
physeq3

BIOMdatamerge_pcoa2 <- ordinate(
  physeq = physeq2, 
  method = "PCoA", 
  distance = "unifrac",  # Change to Unweighted UniFrac
  weighted = FALSE  # Make sure to set 'weighted' to FALSE for Unweighted UniFrac
)


BIOMdatamerge_pcoa2

####unweighted-UniFrac#####
BIOMdatamerge_pcoa3 <- ordinate(
  physeq = physeq2, 
  method = "PCoA", 
  distance = "unweighted-UniFrac",  # Change to Unweighted UniFrac
  weighted = FALSE  # Make sure to set 'weighted' to FALSE for Unweighted UniFrac
)


BIOMdatamerge_pcoa3



# Plot 


p1a = plot_ordination(
  physeq = physeq2,
  ordination = BIOMdatamerge_pcoa2,
  color = "Platform",
  type="samples",
  title = "PCoA of bacterial Communities"
) +
  scale_color_manual(values = c("blue1", "gold", "firebrick1", "black"))+ geom_point(size=3) +
  stat_ellipse()

p1a


###########PlatformandState########
p1a2 <- plot_ordination(
  physeq = physeq2,
  ordination = BIOMdatamerge_pcoa3,
  color = "Year",
  type = "samples",
  title = "PCoA of bacterial Communities"
) +
  theme_bw() +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD")) +
  geom_point(size = 3) +
  scale_fill_manual(values = c("transparent", "black"))

# Add a legend for Platform
p1a2 <- p1a2 + theme(legend.position = "bottom")

# Display the plot
p1a2

p1a4 <- plot_ordination(
  physeq = physeq2,
  ordination = BIOMdatamerge_pcoa3,
  color = "Year",
  type = "samples",
  title = "PCoA of bacterial Communities"
) +
  theme_bw() +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD")) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Platform), type = "norm", level = 0.8, color = "#999999", size = 0.5) +
  scale_fill_manual(values = c("transparent", "black"))

p1a4 <- p1a4 + theme(legend.position = "bottom")
p1a4

p1a5 <- plot_ordination(
  physeq = physeq2,
  ordination = BIOMdatamerge_pcoa3,
  color = "Year",
  type = "samples",
  title = "PCoA of bacterial Communities"
) +
  theme_bw() +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD")) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = State), alpha = 0, size = 1.2) +
  scale_fill_manual(values = c("transparent", "black", "#FFFF99", "#B15928", "#1B9E77"))

p1a5 <- p1a5 + theme(legend.position = "bottom")
p1a5


p1a6 <- plot_ordination(
  physeq = physeq2,
  ordination = BIOMdatamerge_pcoa3,
  color = "Year",
  type = "samples",
  shape= "Platform"
) +
  scale_x_continuous (expand = c(0, 0), breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3), limits = c(-0.4,0.4))+
  theme_bw() +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                "#CAB2D6", "#6A3D9A", "white", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD")) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Platform), type = "norm", level = 0.9, color = "#999999", size = 0.5) + 
  theme(axis.text.x = element_text(size = 14, face='bold'))+
  theme(axis.title=element_text(size=14, face='bold'),
        axis.text=element_text(size=14, face='bold'),
        strip.text=element_text(size=14, face='bold')) + 
  theme(legend.position = "right", 
        plot.title= element_text(size = 14, hjust = 0.5, face='bold')) 

p1a6 <- p1a6 + theme(legend.position = "bottom")
p1a6

plot_richness(data_mergePrelim, x="Year", color="Platform")

plot_richness(data_mergePrelim, measures=c("Chao1", "Shannon"), x="Year", color="Platform") +
  theme_bw() +
  scale_color_manual(values = c( "#FB9A99", "#E31A1C", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FDBF6F", "#FF7F00",
                                 "#CAB2D6", "cadetblue1", "blueviolet","white" ), 
                     breaks = c("Esperance", "Katanning", "Lower_EP", "Upper_EP", "Mid_YP", "Upper_YP", 
                                "Kaniva", "Horsham", "Birchip", "Mid_Kat_Esp", "NSW_outlier")
  ) +
  geom_point(size = 2)


##########
plot_net(physeq2, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.9, color="Phylum", point_label="Class")

###############


#This will
p1a + facet_wrap(~Year, 2) + geom_point(size=3)  + stat_ellipse()

p1a + facet_wrap(~Platform, 2) + geom_point(size=3) + stat_ellipse()



# Plot2
p2 = plot_ordination(physeq = data_merge,
                     ordination = BIOMdatamerge_pcoa,
                     type="samples", color="Year", shape="State") 


p2


p2 + geom_polygon(aes(fill=State)) + geom_point(size=5) + ggtitle("samples")
#LOL

p1a = plot_ordination(
  physeq = data_merge,
  ordination = BIOMdatamerge_pcoa,
  color = "State",
  type="samples",
  title = "PCoA of bacterial Communities"
) +
  scale_color_manual(values = c("blue1", "gold", "firebrick1"))+ geom_point(size=3)


#Include platform

p1a <- plot_ordination(
  physeq = data_merge,
  ordination = BIOMdatamerge_pcoa,
  color = "State",
  type = "samples",
  title = "PCoA of bacterial Communities"
) +
  scale_color_manual(values = c("blue1", "gold", "firebrick1")) +
  geom_point(size = 3) +
  geom_ellipse(aes(fill = Platform), alpha = 0, size = 1.2) +
  scale_fill_manual(values = c("transparent", "black"))  # Adjust ellipse appearance

# Add a legend for Platform
p1a <- p1a + theme(legend.position = "bottom")

# You may need to further customize legend appearance if needed

p1a








p1a + facet_wrap(~Type, 2) + geom_point(size=3)  + stat_ellipse()

p1a + facet_wrap(~State, 2) + geom_point(size=3)


#PCoA by year
# Plot 
p1a2 = plot_ordination(
  physeq = data_merge,
  ordination = BIOMdatamerge_pcoa,
  color = "Year",
  type="samples",
  title = "PCoA of bacterial Communities"
) +
  scale_color_manual(values = c("blue1", "gold", "firebrick1", "black"))+ geom_point(size=3) +
  stat_ellipse()

p1a2

##############Using_rooted_tree#######
Tree2 <- "/Microbiome Stage 1/Prelimstudy_write_up/Export_merged/tree_rooted.nwk" 

# Read the tree file
library("ape")
tree2 <- read.tree(Tree2)

physeq3 = phyloseq(OTU_table_Prelim, psPrelim, mapPrelim, tree2)
physeq3

BIOMdatamerge_pcoa10 <- ordinate(
  physeq = physeq3, 
  method = "PCoA", 
  distance = "unifrac",  # Change to Unweighted UniFrac
  weighted = FALSE  # Make sure to set 'weighted' to FALSE for Unweighted UniFrac
)


BIOMdatamerge_pcoa10

####unweighted-UniFrac#####
BIOMdatamerge_pcoa11 <- ordinate(
  physeq = physeq3, 
  method = "PCoA", 
  distance = "unweighted-UniFrac",  # Change to Unweighted UniFrac
  weighted = FALSE  # Make sure to set 'weighted' to FALSE for Unweighted UniFrac
)


BIOMdatamerge_pcoa11

###PLOTS###

p1a10 <- plot_ordination(
  physeq = physeq3,
  ordination = BIOMdatamerge_pcoa10,
  color = "Year",
  type = "samples",
  title = "PCoA of bacterial Communities"
) +
  theme_bw() +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD")) +
  geom_point(size = 3) +
  scale_fill_manual(values = c("transparent", "black"))

# Add a legend for Platform
p1a10 <- p1a10 + theme(legend.position = "bottom")

# Display the plot
p1a10


p1a11 <- plot_ordination(
  physeq = physeq3,
  ordination = BIOMdatamerge_pcoa10,
  color = "Year",
  type = "samples",
  title = "PCoA of bacterial Communities"
) +
  theme_bw() +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD")) +
  geom_point(size = 3) +
  scale_fill_manual(values = c("transparent", "black"))

# Add a legend for Platform
p1a11 <- p1a11 + theme(legend.position = "bottom")

# Display the plot
p1a11

####HEATMAP#######
plot_heatmap(data_mergePrelimIllumina, method = "PCoA", distance = "bray")

Tree <- "/Microbiome Stage 1/Prelimstudy_write_up/Export_merged/tree.nwk" 

# Read the tree file
library("ape")
tree <- read.tree(Tree)

########Rebuild phyloseq data from scratch using all the simulated data components we just generated###########

physeqIllumina = phyloseq(OTU_table_PrelimIllumina, psPrelimIllumina, mapPrelimIllumina, tax_table(data_mergePrelimIllumina), tree)
physeqIllumina

physeqIllumina_abund <- filter_taxa(physeqIllumina, function(x) sum(x > total*0.01) > 0, TRUE)
physeqIllumina_abund 

PacBio_HM_UniFracIllumina <- plot_heatmap(physeqIllumina_abund, method = "NMDS", distance = "unweighted-UniFrac",
                                        taxa.label = "Class", taxa.order = "Class", low = "lightblue", high = "blue", na.value = "black") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),     # Increase x-axis label font size
    axis.text.y = element_text(size = 12, face = "bold"),     # Increase y-axis label font size
    axis.title = element_text(size = 14, face = "bold"),      # Increase axis title font size
    plot.title = element_text(size = 16, face = "bold")       # Increase plot title font size
  )

PacBio_HM_UniFracIllumina

#"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"

plot_heatmap(data_mergePrelimPacBio, method = "PCoA", distance = "bray")

########Rebuild phyloseq data from scratch using all the simulated data components we just generated###########

physeqPacBio = phyloseq(OTU_table_PrelimPacBio, psPrelimPacBio, mapPrelimPacBio, tax_table(data_mergePrelimPacBio), tree)
physeqPacBio


PacBio_HM_UniFrac <- plot_heatmap(physeqPacBio, method = "PCoA", distance = "unweighted-UniFrac",
             low="lightblue", high="blue", na.value="black")
PacBio_HM_UniFrac

physeqPacBio_abund <- filter_taxa(physeqPacBio, function(x) sum(x > total*0.001) > 0, TRUE)
physeqPacBio_abund


library(ggplot2)

PacBio_HM_UniFracPhylum <- plot_heatmap(physeqPacBio_abund, method = "NMDS", distance = "unweighted-UniFrac",
                                        taxa.label = "Class", taxa.order = "Class", low = "lightblue", high = "blue", na.value = "black", VariableA = c("Year", "State"), 
                                        transformation = "log10", 
                                        annotation_colors=meta_colors) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),     # Increase x-axis label font size
    axis.text.y = element_text(size = 12, face = "bold"),     # Increase y-axis label font size
    axis.title = element_text(size = 14, face = "bold"),      # Increase axis title font size
    plot.title = element_text(size = 16, face = "bold")       # Increase plot title font size
  )

PacBio_HM_UniFracPhylum

##########FUN with Pheatmap#########
library(microbiome)
library(microbiomeutilities)
library(knitr)
library(tibble)
library(dplyr)
library(viridis)
library(reshape2)
#> Loading required package: viridisLite
library(RColorBrewer)

kable(head(tax_table(physeqPacBio)))

physeqPacBio2 <- core(physeqPacBio, detection = 10, prevalence = 20 / 100)
kable(head(tax_table(physeqPacBio2)))


meta_colors <- list(c("2001" = "#A6CEE3", 
                      "2002" = "#1F78B4",
                      "2003" = "#B2DF8A", 
                      "2004" = "#33A02C",
                      "2006" = "#FB9A99", 
                      "2008" = "#E31A1C",
                      "2011" = "#FDBF6F", 
                      "2014" = "#FF7F00",
                      "2017" = "#CAB2D6", 
                      "2021" = "#6A3D9A"
                      ), 
                    c("VIC" = "darkseagreen1", 
                      "SA" = "#FFFF99", 
                      "WA"="lightpink"))


# add labels for pheatmap to detect
names(meta_colors) <- c("Year", "State")

# create a gradient color palette for abundance
grad_ab <- colorRampPalette(c("#faf3dd","#f7d486" ,"#5e6472"))
grad_ab_pal <- grad_ab(10)

physeqPacBio_abund <- filter_taxa(physeqPacBio, function(x) sum(x > total*0.001) > 0, TRUE)


heat.sample3 <- plot_taxa_heatmap(physeqPacBio_abund,
                                 subset.top = 60,
                                 VariableA = c("Year", "State"), 
                                 heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                 transformation = "log10", 
                                 cluster_rows = T, cluster_coloums = F,
                                 show_colnames = T,
                                 annotation_colors=meta_colors) 

heat.sample3


heat.sample4 <- plot_taxa_heatmap(physeqIllumina_abund,
                                  subset.top = 60,
                                  VariableA = c("Year", "State"), 
                                  heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                  transformation = "log10", 
                                  cluster_rows = T, cluster_coloums = F,
                                  show_colnames = T,
                                  annotation_colors=meta_colors) 

heat.sample4
data_mergePrelim



# Remove samples with specific names
samples_to_remove <- c("NTC_2.fastq", "ZymoBIOMICS_Microbial_Community_DNA_Standard_1.fastq", "ZymoBIOMICS_Microbial_Community_DNA_Standard_2.fastq", 
                       "ZymoBIOMICS_Microbial_Community_DNA_Standard_3.fastq", "ZymoBIOMICS_Microbial_Community_DNA_Standard_4.fastq")
data_mergePrelim_filtered <- subset_samples(physeq3, !(sample_names(physeq3) %in% samples_to_remove))

data_mergePrelim_filtered

#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(data_mergePrelim_filtered))
standf = function(x, t=total) round(t * (x / sum(x)))
data_mergePrelim_filtered = transform_sample_counts(data_mergePrelim_filtered, standf)
data_mergePrelim_filtered

data_mergePrelimabund <- filter_taxa(data_mergePrelim_filtered, function(x) sum(x > total*0.015) > 0, TRUE)
data_mergePrelimabund



plot_heatmap(data_mergePrelimabund, method = "PCoA", distance = "unweighted-UniFrac")



# Generate the heatmap
plot_heatmap(data_mergePrelim_filtered, method = "PCoA", distance = "unweighted-UniFrac", 
             low="lightblue", high="blue", na.value="black")


samples_to_remove <- c("2023_PacBio")
data_mergePrelim_filtered <- subset_samples(physeq3, !(sample_variables(physeq3) %in% samples_to_remove))

 
plot_richness(data_mergePrelim2, measures=c("Observed", "Chao1", "ACE", "Simpson"), x="Platform", color="Year") +
  theme_bw() +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                 "#CAB2D6", "blueviolet","cadetblue1", "white" )
  ) +
  geom_point(size = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 13, face = 'bold'),
        strip.text = element_text(size = 13, face = 'bold'),
        strip.background = element_rect(fill = NA, color = 'black', linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = 'black'),
        legend.position = 'right',
        plot.title = element_text(size = 18, hjust = 0.5))


# Extract richness data
richness_data <- estimate_richness(data_mergePrelim_filtered, measures = c("Observed", "Chao1", "ACE", "Simpson", "invsimpson"))
head(richness_data)

# Convert sample data to data frame
sample_data_df <- as.data.frame(sample_data(data_mergePrelim_filtered))
head(sample_data_df)


# Convert richness data to a data frame
richness_df <- as.data.frame(richness_data)

Illumina_PacBio_Alpha_data_ratios <- read_excel("Desktop/Microbiome Stage 1/Prelimstudy_write_up/Illumina_PacBio_merge_for_anal/Illumina_PacBio_Alpha_data_ratios.xlsx")

 
data <- Illumina_PacBio_Alpha_data_ratios

library(gridExtra)

# Create individual box plots
boxplot_observed <- ggplot(data, aes(x = as.factor(Year), y = Observed, color = as.factor(Year))) + 
  geom_boxplot() + 
  scale_color_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                              "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                              "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD"))+
  labs(title = "Observed", x="Year", y = "Ratio Observed", color = "Year") +
  scale_y_continuous(expand = c(0, 0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8), limits = c(0,1.8)) +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1, face='bold'),
    axis.title = element_text(size=12, face='bold'),
    axis.text = element_text(size=12, face='bold'),
    strip.text = element_text(size=12, face='bold'),
    strip.background = element_rect(fill=NA, color='black', size=0.5),
    panel.border = element_rect(fill=NA, color='black'),
    legend.position = "right", 
    plot.title = element_text(size = 14, face='bold', hjust = 0.5, vjust = 0.5)
  ) 

boxplot_chao1 <- ggplot(data, aes(x = as.factor(Year), y = Chao1, color = as.factor(Year))) + 
  geom_boxplot() + 
  labs(title = "Chao1", x="Year", y = "Ratio Chao1", color = "Year") +
  scale_color_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                              "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                              "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD"))+
  scale_y_continuous(expand = c(0, 0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), limits = c(0,2.0)) +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1, face='bold'),
    axis.title = element_text(size=12, face='bold'),
    axis.text = element_text(size=12, face='bold'),
    strip.text = element_text(size=12, face='bold'),
    strip.background = element_rect(fill=NA, color='black', size=0.5),
    panel.border = element_rect(fill=NA, color='black'),
    legend.position = "right", 
    plot.title = element_text(size = 14, face='bold', hjust = 0.5, vjust = 0.5)
  ) 

boxplot_ace <- ggplot(data, aes(x = as.factor(Year), y = ACE, color = as.factor(Year))) + 
  geom_boxplot() + 
  scale_color_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                              "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                              "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD"))+
  labs(title = "ACE", x="Year", y = "Ratio ACE", color = "Year") +
  scale_y_continuous(expand = c(0, 0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0), limits = c(0,2.0)) +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1, face='bold'),
    axis.title = element_text(size=12, face='bold'),
    axis.text = element_text(size=12, face='bold'),
    strip.text = element_text(size=12, face='bold'),
    strip.background = element_rect(fill=NA, color='black', size=0.5),
    panel.border = element_rect(fill=NA, color='black'),
    legend.position = "right", 
    plot.title = element_text(size = 14, face='bold', hjust = 0.5, vjust = 0.5)
  ) 

boxplot_simpson <- ggplot(data, aes(x = as.factor(Year), y = Simpson, color = as.factor(Year))) + 
  geom_boxplot() + 
  scale_color_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                              "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                              "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD"))+
  labs(title = "Simpson", x="Year", y = "Ratio Simpson", color = "Year") +
  scale_y_continuous(expand = c(0, 0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), limits = c(0,1.2)) +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1, face='bold'),
    axis.title = element_text(size=12, face='bold'),
    axis.text = element_text(size=12, face='bold'),
    strip.text = element_text(size=12, face='bold'),
    strip.background = element_rect(fill=NA, color='black', size=0.5),
    panel.border = element_rect(fill=NA, color='black'),
    legend.position = "right", 
    plot.title = element_text(size = 14, face='bold', hjust = 0.5, vjust = 0.5)
  ) 

# Combine them in a 2x2 grid
boxplots_grid <- grid.arrange(boxplot_observed, boxplot_chao1, boxplot_ace, boxplot_simpson, ncol = 2)


########Hill diversity#####
# Calculate Simpson index
simpson_index <- function(abundance_vector) {
  N <- sum(abundance_vector)
  p_i <- abundance_vector / N
  return(1 - sum(p_i^2))
}

# Apply the function to each sample
simpson_values <- apply(otu_table(data_mergePrelim_filtered), 2, simpson_index)

# Create a data frame with sample IDs and Simpson indices
simpson_df <- data.frame(Sample = colnames(otu_table(data_mergePrelim_filtered)), Simpson_Index = simpson_values)
simpson_df

# Create a bar plot of Simpson indices
ggplot(simpson_df, aes(x = Sample, y = Simpson_Index)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(x = "Sample", y = "Simpson Index") +
  theme_minimal()

boxplot_simpson_Hill <- ggplot(simpson_df, aes(x = Sample, y = Simpson_Index)) + 
  geom_boxplot() + 
  scale_color_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                              "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                              "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD"))+
  labs(title = "Simpson", x="Year", y = "Ratio Simpson", color = "Year") +
  scale_y_continuous(expand = c(0, 0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), limits = c(0,1.2)) +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1, face='bold'),
    axis.title = element_text(size=12, face='bold'),
    axis.text = element_text(size=12, face='bold'),
    strip.text = element_text(size=12, face='bold'),
    strip.background = element_rect(fill=NA, color='black', size=0.5),
    panel.border = element_rect(fill=NA, color='black'),
    legend.position = "right", 
    plot.title = element_text(size = 14, face='bold', hjust = 0.5, vjust = 0.5)
  ) 

boxplot_simpson_Hill

########Diversity 16S#########

#' Diversity function for 16S amplicon data
#'
#' This function calculates Hill diversity metrics from 16S amplicon data (in phyloseq format).
#' D0 (richness) is calculated with three methods: 1) Observed richness, 2) Chao1 estimator
#' 3) breakaway method (Willis and Bunge 2016). D1 (exponential of Shannon entropy)
#' and D2 (inverse Simpson index) are respectively Hill order 1 and 2.
#' Errors for D1 and D2 are calculated by bootstrapping.
#' 
#' @param x phyloseq object generated by the phyloseq package
#' @param R Number of bootstraps to conduct. Defaults to 999
#' @param brea TRUE/FALSE if breakaway method for D0 estimation should be used.
#' Defaults to TRUE. This method fails easily if you don't have atleast 6 contiguous
#' frequencies.
#' @param thresh Minimum sample size required to perform Chao1 estimation.
#' @param parallel Should the calculation be parallelized? Defaults to FALSE
#' @param ncores How many cores should be used in case of parallel computation?
#' Defaults to 2.
#' @keywords diversity, fcm, alpha
#' @importFrom foreach %dopar%
#' @importFrom stats sd
#' @examples 
#' # First install phyloseq if you haven't yet
#' if(requireNamespace("phyloseq",quietly=TRUE)){
#' library(phyloseq)
#' } else {
#' source("https://bioconductor.org/biocLite.R")
#' biocLite("phyloseq")
#' library(phyloseq)
#' }
#' # Load data (V3-V4 amplicon data from doi: 10.1111/2041-210X.12607)
#' data(physeq_test)
#' # Opting for three bootstraps, because this can take some time.
#' Diversity_16S(phyloseq::prune_samples(phyloseq::sample_names(physeq_test) == "1", physeq_test), R=3)
#' @export

Diversity_16S <- function(x, R = 999, brea = TRUE, thresh = 200, parallel = FALSE, 
                          ncores = 2) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Phyloseq package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  cat("\t**WARNING** this functions assumes that rows are samples and columns
      \tare taxa in your phyloseq object, please verify.\n")
  
  # Matrix for storing data
  DIV <- matrix(nrow = phyloseq::nsamples(x), ncol = 10)
  row.names(DIV) <- phyloseq::sample_names(x)
  
  # Diversity functions
  D0.boot <- function(x) sum(x != 0)
  D2.boot <- function(x) 1/sum((x)^2)
  D1.boot <- function(x) exp(-sum(x * log(x)))
  
  if(parallel == FALSE){
    # Start resampling
    for (i in 1:phyloseq::nsamples(x)) {
      temp.D0 <- c()
      temp.D1 <- c()
      temp.D2 <- c()
      temp.phy <- phyloseq::prune_samples(x = x, samples = phyloseq::sample_names(x)[i])
      cat(paste0(date(), "\tCalculating diversity for sample ",i,"/",phyloseq::nsamples(x)," --- ",phyloseq::sample_names(x)[i], "\n"))
      for (j in 1:R) {
        temp <- phyloseq::rarefy_even_depth(temp.phy, verbose = FALSE, replace = TRUE)
        # Calculate frequencies
        temp <- data.frame(phyloseq::transform_sample_counts(temp, fun = function(x) x/sum(x))@otu_table)
        # Calculate Diversities
        temp.D0 <- c(temp.D0, D0.boot(temp))
        temp.D1 <- c(temp.D1, D1.boot(temp))
        temp.D2 <- c(temp.D2, D2.boot(temp))
        # Store diversities at the end of resampling run
        if (j == R) {
          DIV[i, 1] <- mean(temp.D0)
          DIV[i, 2] <- stats::sd(temp.D0)
          DIV[i, 7] <- mean(temp.D1)
          DIV[i, 8] <- stats::sd(temp.D1)
          DIV[i, 9] <- mean(temp.D2)
          DIV[i, 10] <- stats::sd(temp.D2)
          remove(temp.D0, temp.D1, temp.D2)
          # cat(paste0(date(), "\tDone with sample ", phyloseq::sample_names(x)[i], "\n"))
        }
      }
      # Perform breakaway for richness estimation
      temp <- t(matrix(temp.phy@otu_table))
      temp <- temp[temp != 0]
      temp <- data.frame(table(temp))
      temp <- apply(temp, 2, FUN = function(x) as.integer(x))
      if(brea==TRUE && phyloseq::sample_sums(temp.phy) > thresh){
        rich <- breakaway::breakaway(temp, print = FALSE, plot = FALSE, answers = TRUE, force=TRUE)
        if(!is.null(rich)){
          DIV[i, 3] <- rich$est
          DIV[i, 4] <- rich$seest
        } else {
          DIV[i, 3] <- NA
          DIV[i, 4] <- NA
        }
      }
      if(phyloseq::sample_sums(temp.phy) >= thresh){
        rich.chao <- breakaway::chao1(temp, output = FALSE, answers = TRUE)
        DIV[i, 5] <- rich.chao$est
        DIV[i, 6] <- rich.chao$seest
      } else {
        DIV[i, 5] <- NA
        DIV[i, 6] <- NA
      }
    }
  } else {
    # Initiate/register multiple cores
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    cat(date(), "\tUsing", ncores, "cores for calculations\n")
    # Start resampling
    for (i in 1:phyloseq::nsamples(x)) {
      temp.D0 <- c()
      temp.D1 <- c()
      temp.D2 <- c()
      temp.phy <- phyloseq::prune_samples(x = x, samples = phyloseq::sample_names(x)[i])
      cat(paste0(date(), "\tCalculating diversity for sample ", i, "/", phyloseq::nsamples(x)," --- ", phyloseq::sample_names(x)[i], "\n"))
      
      # Parallelize diversity calculations 
      tmp <- foreach::foreach(j = 1:R, .combine = rbind) %dopar% {
        temp <- phyloseq::rarefy_even_depth(temp.phy, verbose = FALSE, replace = TRUE)
        # Calculate frequencies
        temp <- data.frame(phyloseq::transform_sample_counts(temp, fun = function(x) x/sum(x))@otu_table)
        # Calculate Diversities
        cbind(D0.boot(temp), D1.boot(temp), D2.boot(temp))
      }
      
      DIV[i, c(1,7,9)] <- colMeans(tmp)
      DIV[i, c(2,8,10)] <- apply(tmp, 2, stats::sd)
      # if (i > 1) DIV <- rbind(DIV, matrix(c(colMeans(tmp), apply(tmp, 2, FUN = sd)), nrow = 1))
      
      # Perform breakaway for richness estimation
      temp <- t(matrix(temp.phy@otu_table))
      temp <- temp[temp != 0]
      temp <- data.frame(table(temp))
      temp <- apply(temp, 2, FUN = function(x) as.integer(x))
      if(brea==TRUE && phyloseq::sample_sums(temp.phy) > thresh){
        rich <- breakaway::breakaway(temp, print = FALSE, plot = FALSE, answers = TRUE, force=TRUE)
        if(!is.null(rich)){
          DIV[i, 3] <- rich$est
          DIV[i, 4] <- rich$seest
        } else {
          DIV[i, 3] <- NA
          DIV[i, 4] <- NA
        }
      }
      if(phyloseq::sample_sums(temp.phy) >= thresh){
        rich.chao <- breakaway::chao1(temp, output = FALSE, answers = TRUE)
        DIV[i, 5] <- rich.chao$est
        DIV[i, 6] <- rich.chao$seest
      } else {
        DIV[i, 5] <- NA
        DIV[i, 6] <- NA
      }
      cat(paste0(date(), "\tDone with sample ", phyloseq::sample_names(x)[i], "\n"))
    }
    
    if(parallel == TRUE){
      cat(date(), "\tClosing workers\n")
      parallel::stopCluster(cl)
    }
  }
  
  colnames(DIV) = c("D0", "sd.D0", "D0.bre" , "sd.D0.bre", "D0.chao", "sd.D0.chao", "D1", "sd.D1", "D2", 
                    "sd.D2")
  cat(date(), "\tDone with all", phyloseq::nsamples(x), "samples\n")
  return(DIV)
}

library(breakaway)

Hill_diversity_999perm <- Diversity_16S(
  data_mergePrelim_filtered,
  R = 999,
  brea = TRUE,
  thresh = 200,
  parallel = FALSE,
  ncores = 2
)


#D0 (Mean of non-zero count data)
#sd.D0 (Standard deviation of non-zero count data)
#D0.bre (Breakaway richness estimate)
#sd.D0.bre (Standard deviation of Breakaway richness estimate)
#D0.chao (Chao1 richness estimate)
#sd.D0.chao (Standard deviation of Chao1 richness estimate)
#D1 (Shannon diversity index)
#sd.D1 (Standard deviation of Shannon diversity index)
#D2 (Inverse Simpson diversity index)
#sd.D2 (Standard deviation of Inverse Simpson diversity index)

#########Manual Hill with coverage#########
# Define the Shannon index function
shannon_index <- function(abundance_vector) {
  p_i <- abundance_vector / sum(abundance_vector)
  p_i[p_i == 0] <- NA  # Exclude zero abundances from the calculation
  return(-sum(p_i * log(p_i), na.rm = TRUE))
}

# Apply the function to each sample
shannon_values <- apply(otu_table(data_mergePrelim_filtered), 2, shannon_index)

# Create a data frame with sample IDs and Shannon indices
shannon_df <- data.frame(Sample = colnames(otu_table(data_mergePrelim_filtered)), Shannon_Index = shannon_values)
shannon_df


##########Shannon HILL With_filered_for chlorplasts_etc_and ngeg_controls removed########
phy.Prelim

# Apply the function to each sample
shannon_values <- apply(otu_table(phy.Prelim), 2, shannon_index)

# Create a data frame with sample IDs and Shannon indices
shannon_df <- data.frame(Sample = colnames(otu_table(phy.Prelim)), Shannon_Index = shannon_values)
shannon_df

# export data from df
write.csv(shannon_df, file = "/Prelim study/Extra figures/df_melt_unrarified_data_ShannonHill.csv", row.names = FALSE)


#Rarefied_data
r1.16s_Prelim

# Apply the function to each sample
shannon_valuesR <- apply(otu_table(r1.16s_Prelim), 2, shannon_index)

# Create a data frame with sample IDs and Shannon indices
shannon_dfrare <- data.frame(Sample = colnames(otu_table(r1.16s_Prelim)), Shannon_Index = shannon_valuesR)
shannon_dfrare

write.csv(shannon_dfrare, file = "/Prelim study/Extra figures/df_melt_rarified_data_ShannonHill.csv", row.names = FALSE)


##############NEW_Phyloseq_objects_for_illumina_PacBio############
########Illumina##########

# Remove samples with specific names
# Set the path to the BIOM file
biom_file <- "/Microbiome Stage 1/Prelimstudy_write_up/Export_merged/table-with-taxonomy.biom"

# Import the BIOM file into a phyloseq object
psPrelimIllumina <- import_biom(biom_file)
psPrelimIllumina

#Check tax labels in BIOM
colnames(tax_table(psPrelimIllumina))

#Check sample labels in BIOM
colnames(otu_table(psPrelimIllumina))

#rename tax labels
colnames(tax_table(psPrelimIllumina)) <- c("Kingdom", "Phylum", "Class", 
                                    "Order", "Family", "Genus", "Species")

# Import sample metadata
mapPrelimIllumina <- read_xlsx("/Microbiome Stage 1/Prelimstudy_write_up/Manifest_All_prelim_study_forfeturetable_testIllumina.xlsx")
mapPrelimIllumina


#convert year to a factor
mapPrelimIllumina[, 'Year'] <- lapply(mapPrelimIllumina[, 'Year'], factor)
mapPrelimIllumina

mapPrelimIllumina <- sample_data(mapPrelimIllumina)
mapPrelimIllumina


#Assign row names to be Sample ID's
row.names(mapPrelimIllumina) <- mapPrelimIllumina$id
mapPrelimIllumina

sample_names(mapPrelimIllumina)

sample_names(psPrelimIllumina)


# Merge data object with sample metadata
data_mergePrelimIllumina <- merge_phyloseq(psPrelimIllumina, mapPrelimIllumina)
data_mergePrelimIllumina

sample_names(data_mergePrelimIllumina)

colnames(tax_table(data_mergePrelimIllumina))


OTU_table_PrelimIllumina <- otu_table(data_mergePrelimIllumina)


############PacBio############

# Set the path to the BIOM file
biom_file <- "/Microbiome Stage 1/Prelimstudy_write_up/Export_merged/table-with-taxonomy.biom"

# Import the BIOM file into a phyloseq object
psPrelimPacBio <- import_biom(biom_file)
psPrelimPacBio

#Check tax labels in BIOM
colnames(tax_table(psPrelimPacBio))

#Check sample labels in BIOM
colnames(otu_table(psPrelimPacBio))

#rename tax labels
colnames(tax_table(psPrelimPacBio)) <- c("Kingdom", "Phylum", "Class", 
                                    "Order", "Family", "Genus", "Species")

# Import sample metadata
mapPrelimPacBio <- read_xlsx("/Microbiome Stage 1/Prelimstudy_write_up/Manifest_All_prelim_study_forfeturetable_testPacBio.xlsx")
mapPrelimPacBio


#convert year to a factor
mapPrelimPacBio[, 'Year'] <- lapply(mapPrelimPacBio[, 'Year'], factor)
mapPrelimPacBio

mapPrelimPacBio <- sample_data(mapPrelimPacBio)
mapPrelimPacBio


#Assign row names to be Sample ID's
row.names(mapPrelimPacBio) <- mapPrelimPacBio$id
mapPrelimPacBio

sample_names(mapPrelimPacBio)

sample_names(psPrelimPacBio)


# Merge data object with sample metadata
data_mergePrelimPacBio <- merge_phyloseq(psPrelimPacBio, mapPrelimPacBio)
data_mergePrelimPacBio

sample_names(data_mergePrelimPacBio)

colnames(tax_table(data_mergePrelimPacBio))


OTU_table_PrelimPacBio <- otu_table(data_mergePrelimPacBio)
head(OTU_table_PrelimPacBio)


ps_venn(
  data_mergePrelimIllumina,
  "State", fraction = 0.5, quantities = list(type=c("percent", "counts")), relative = TRUE, font = 2, 
  weight = FALSE, labels = list(cex = 2), col = "black", fill = c("pink","blue","purple", "orange"))


ps_venn(
  data_mergePrelimPacBio,
  "State", fraction = 0.5, quantities = list(type=c("percent", "counts")), relative = TRUE, font = 2, 
  weight = FALSE, labels = list(cex = 2), col = "black", fill = c("pink","blue","purple", "orange"))

ps_venn(
  data_mergePrelim,
  "Platform", fraction = 0.001, quantities = list(type=c("percent", "counts")), relative = TRUE, font = 2, 
  weight = FALSE, labels = list(cex = 2), col = "black", fill = c("pink","blue","purple", "orange"))


#######################################################
############Filtered_alpha_rare_calcs###################

#Remove negative controls
phy_inStates <- subset_samples(data_mergePrelim, State == "VIC" | State == "WA" | State == "SA")
phy_inStates

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 43460 taxa and 173 samples ]
#sample_data() Sample Data:       [ 173 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 43460 taxa by 7 taxonomic ranks ]

min( sample_sums(phy_inStates) ) # 85
min(taxa_sums(phy_inStates)) # 0

# prune taxa that have zero sequence reads
phy_inStates <- prune_taxa(taxa = taxa_sums(phy_inStates) > 0, x = phy_inStates)
phy_inStates
# phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 43372 taxa and 173 samples ]
#sample_data() Sample Data:       [ 173 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 43372 taxa by 7 taxonomic ranks ]

min(sample_sums(phy_inStates)) # 85
min(taxa_sums(phy_inStates)) # 1


#Filter out eukaryotes and mitochondria
#First Check the format taxa names are in
rank_names(phy_inStates)
tax_table(phy_inStates)

######use this for filtering!!########
#Order instead of class for chloroplast!#########

phy.Prelim <- phy_inStates %>%
  subset_taxa(
    Kingdom == "d__Bacteria" &                   #only bacteria
      Family  != "f__Mitochondria" &             #filter out mitochondria
      Order   != "o__Chloroplast"                #filter out chloroplasts
  )
phy.Prelim

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 38572 taxa and 173 samples ]
#sample_data() Sample Data:       [ 173 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 38572 taxa by 7 taxonomic ranks ]


################
# Load the DT library
library(DT)
library(picante)
library(microbiome) # data analysis and visualisation
library(microbiomeutilities) # some utility tools 
library(ggpubr) # publication quality figures, based on ggplot2
library(data.table) # alternative to data.frame
library(dplyr) # data handling  

sdt = data.table(as(sample_data(phy.Prelim), "data.frame"),
                 TotalReads = sample_sums(phy.Prelim), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")

pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram(color = "black", fill = "#7AA6DCFF", binwidth = 500) + ggtitle("Sequencing Depth") +
  theme_bw() 
pSeqDepth 

pSeqDepth + facet_wrap(~Year)

library(ggplot2)

ggplot(sdt, aes(Year, TotalReads, color = Year, alpha = Year)) + 
  geom_point(size = 2) +
  geom_jitter(size = 2, width = 0.5) +
  ggtitle("Sequencing Depth vs. Year") +
  scale_y_log10() +
  theme_bw() +
  scale_color_manual(values = c("#A73030FF", "#EFC000FF", "#0073C2FF", "#868686FF", "#8F7700FF", "#3B3B3BFF", "#7876B1FF", "#FFDC91FF", "#E64B35FF", "#4DBBD5FF",
                                "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#20854EFF", "#E18727FF", 
                                "#B1746FFF", "#ADB17DFF",   "#EE4C97FF","#7AA6DCFF")) +
  scale_alpha_manual(values = c(YrRecieved = 0.05)) +  # Set alpha values for each level
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  guides(alpha = guide_legend(title = "Year Received")) +  # Customize the legend title
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  facet_wrap(~Platform)

sdt[, min(TotalReads)] #83

########remove low read counts#########
sample_sums(phy.Prelim)

summary(sample_sums(phy.Prelim))

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#83    3490   52271   44488   83983  125832 

#want to keep samples with at least 1500 sequence reads
min_reads_threshold <- 1500

# Filter samples
phy.Prelim_filtered <- prune_samples(sample_sums(phy.Prelim) >= min_reads_threshold, phy.Prelim)
phy.Prelim_filtered 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 38572 taxa and 172 samples ]
#sample_data() Sample Data:       [ 172 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 38572 taxa by 7 taxonomic ranks ]

##lost one sample AAE0543

summary(sample_sums(phy.Prelim_filtered  ))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1563    3492   52478   44746   84247  125832  

sdt = data.table(as(sample_data(phy.Prelim_filtered ), "data.frame"),
                 TotalReads = sample_sums(phy.Prelim_filtered ), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")

ggplot(sdt, aes(Year, TotalReads, color = Year, alpha = Year)) + 
  geom_point(size = 2) +
  geom_jitter(size = 2, width = 0.5) +
  ggtitle("Sequencing Depth vs. Year") +
  scale_y_log10() +
  theme_bw() +
  scale_color_manual(values = c("#A73030FF", "#EFC000FF", "#0073C2FF", "#868686FF", "#8F7700FF", "#3B3B3BFF", "#7876B1FF", "#FFDC91FF", "#E64B35FF", "#4DBBD5FF",
                                "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#20854EFF", "#E18727FF", 
                                "#B1746FFF", "#ADB17DFF",   "#EE4C97FF","#7AA6DCFF")) +
  scale_alpha_manual(values = c(YrRecieved = 0.05)) +  # Set alpha values for each level
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  guides(alpha = guide_legend(title = "Year Received")) +  # Customize the legend title
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  facet_wrap(~Platform)


########Checking original data post trim##########
##need to seperate PacBio from Illumina
phy_inPacBio <- subset_samples(phy.Prelim_filtered, Platform == "PacBio")
phy_inPacBio
phy_inIllumina <- subset_samples(phy.Prelim_filtered, Platform == "Illumina")
phy_inIllumina

lib.divIll <- alpha(phy_inPacBio, index = "all")
lib.divIll <- alpha(phy_inIllumina, index = "all")



#PacBio
# let us add nmber of total reads/samples
lib.div2 <- richness(phy_inPacBio)

lib.divIll$ReadsPerSample <- sample_sums(phy_inPacBio)

lib.divIll$Observed <- lib.div2$observed

colnames(lib.divIll)

#Use ggscatter function from ggpubr package to visualze.

p1 <- ggscatter(lib.divIll, "diversity_shannon", "ReadsPerSample") +
  stat_cor(method = "pearson")

p2 <- ggscatter(lib.divIll, "diversity_inverse_simpson", "ReadsPerSample",
                add = "loess"
) +
  stat_cor(method = "pearson")

p3 <- ggscatter(lib.divIll, "Observed", "ReadsPerSample",
                add = "loess") +
  stat_cor(
    method = "pearson",
    label.x = 100,
    label.y = 50000
  )

ggarrange(p1, p2, p3, ncol = 2, nrow = 2)

#Illumina
# let us add nmber of total reads/samples
lib.div2 <- richness(phy_inIllumina)

lib.divIll$ReadsPerSample <- sample_sums(phy_inIllumina)

lib.divIll$Observed <- lib.div2$observed

colnames(lib.divIll)

#Use ggscatter function from ggpubr package to visualze.

p1 <- ggscatter(lib.divIll, "diversity_shannon", "ReadsPerSample") +
  stat_cor(method = "pearson")

p2 <- ggscatter(lib.divIll, "diversity_inverse_simpson", "ReadsPerSample",
                add = "loess"
) +
  stat_cor(method = "pearson")

p3 <- ggscatter(lib.divIll, "Observed", "ReadsPerSample",
                add = "loess") +
  stat_cor(
    method = "pearson",
    label.x = 100,
    label.y = 50000
  )

ggarrange(p1, p2, p3, ncol = 2, nrow = 2)



#######Checking rarefied data########
####rarefy with higher seq count########
seed <- 1234
r1.16s_Prelim <- rarefy_even_depth(phy.Prelim_filtered, sample.size = min(sample_sums(phy.Prelim_filtered)),
                                    rngseed = seed, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

min(taxa_sums(r1.16s_Prelim)) # 1
sample_sums(r1.16s_Prelim) # all 1563
ntaxa(r1.16s_Prelim) # 30201

#########check taxa prevalence###############

r1.16s_Prelim_taxa <- plot_taxa_prevalence(r1.16s_Prelim, "Phylum") +
  theme(legend.text = element_text(size = 12)) 

r1.16s_Prelim_taxa


####Rarified seq depth###########
sdt = data.table(as(sample_data(r1.16s_Prelim), "data.frame"),
                 TotalReads = sample_sums(r1.16s_Prelim), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")

ggplot(sdt, aes(Year, TotalReads, color = Year, alpha = Year)) + 
  geom_point(size = 2) +
  geom_jitter(size = 2, width = 0.5) +
  ggtitle("Sequencing Depth vs. Time") +
  scale_y_log10() +
  theme_bw() +
  scale_color_manual(values = c("#A73030FF", "#EFC000FF", "#0073C2FF", "#868686FF", "#8F7700FF", "#3B3B3BFF", "#7876B1FF", "#FFDC91FF", "lightpink", "#4DBBD5FF",
                                "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#7E6148FF", "thistle", "#E18727FF", 
                                "#B1746FFF", "slateblue",   "palevioletred","#7AA6DCFF")) +
  scale_alpha_manual(values = c(YrRecieved = 0.05)) +  # Set alpha values for each level
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  guides(alpha = guide_legend(title = "Year Received")) +  # Customize the legend title
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  facet_wrap(~Platform)


#####Shannon Diversity Index#########
shan <- plot_richness(r1.16s_Prelim, measures=c("Shannon"))

out <- data.frame(sample=shan$data$id,
                  shannon=shan$data$value,
                  Year=shan$data$Year,
                  Platform=shan$data$Platform)

out$Eff_No_ASV <- exp(out$shannon)

names(out)
# [1] "sample"     "shannon"    "Year"       "Platform"   "Eff_No_ASV"

library(data.table)
df.melt <- melt(as.data.table(out[, c("sample", "Year", "Platform", "Eff_No_ASV")]), id.vars = c("sample", "Year", "Platform"))
df.melt

str(df.melt)
# Classes ‘data.table’ and 'data.frame':	172 obs. of  4 variables:
#$ sample  : chr  "AAA5173.fastq" "AAA5228.fastq" "AAA5670.fastq" "AAA5727.fastq" ...
#$ Year    : Factor w/ 10 levels "2001","2002",..: 1 1 1 1 1 1 1 1 3 3 ...
#$ Platform: chr  "PacBio" "PacBio" "PacBio" "PacBio" ...
#$ variable: Factor w/ 1 level "Eff_No_ASV": 1 1 1 1 1 1 1 1 1 1 ...
#$ value   : num  498 371 378 231 434 ...
#- attr(*, ".internal.selfref")=<externalptr>   


unique(df.melt$Year)
# [1] 2001 2003 2002 2004 2006 2008 2011 2014 2017 2021
#Levels: 2001 2002 2003 2004 2006 2008 2011 2014 2017 2021 

df.melt$Year <- factor(df.melt$Year,
                       levels=c(
                         "2001", "2002", "2003", "2004", "2006", 
                         "2008", "2011", "2014", "2017", "2021"
                       ),
                       labels=c(
                         "2001", "2002", "2003", "2004", "2006", 
                         "2008", "2011", "2014", "2017", "2021"
                       ),
                       ordered = TRUE)

cols <- c("2001" = "#A73030FF", 
          "2002" = "#EFC000FF",
          "2003" = "#0073C2FF",
          "2004" = "#868686FF",
          "2006" = "#3B3B3BFF",
          "2008" = "#7876B1FF",
          "2011" = "#4DBBD5FF",
          "2014" = "#F39B7FFF",
          "2017" =  "#7E6148FF",
          "2021" =  "slateblue"
)

pr1.16s_Prelim <- ggplot(data=df.melt, aes(x=Year, y=value, color = Year)) +
  ggtitle("Bacteria alpha diversity ASV's")+
  theme_bw() +
  #geom_point() +
  geom_boxplot() +
  geom_jitter(size=1.5, width = 0.15) +
  scale_colour_manual(values=cols)+
  theme(axis.text.x  = element_text(angle=90, hjust=1, vjust = 1) ) +
  #facet_wrap("variable",scales = "free", nrow = 1, ncol = 3) +
  labs(x = NULL, y = "Effective no of ASVs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  guides(alpha = guide_legend(title = "Year Received")) +  # Customize the legend title
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  facet_wrap(~Platform)
pr1.16s_Prelim


df.melt2 <- melt(as.data.table(out[, c("sample", "Year", "Platform", "shannon")]), id.vars = c("sample", "Year", "Platform"))

str(df.melt2)
# Classes ‘data.table’ and 'data.frame':	172 obs. of  5 variables:
#$ sample  : chr  "AAA5173.fastq" "AAA5228.fastq" "AAA5670.fastq" "AAA5727.fastq" ...
#$ Year    : Factor w/ 10 levels "2001","2002",..: 1 1 1 1 1 1 1 1 3 3 ...
#$ Platform: chr  "PacBio" "PacBio" "PacBio" "PacBio" ...
#$ variable: Factor w/ 1 level "shannon": 1 1 1 1 1 1 1 1 1 1 ...
#$ value   : num  6.21 5.92 5.94 5.44 6.07 ...
#- attr(*, ".internal.selfref")=<externalptr>  

df.melt2$Year <- factor(df.melt2$Year,
                       levels=c(
                         "2001", "2002", "2003", "2004", "2006", 
                         "2008", "2011", "2014", "2017", "2021"
                       ),
                       labels=c(
                         "2001", "2002", "2003", "2004", "2006", 
                         "2008", "2011", "2014", "2017", "2021"
                       ),
                       ordered = TRUE)

cols <- c("2001" = "#A6CEE3", 
          "2002" = "#1F78B4",
          "2003" = "#B2DF8A",
          "2004" = "#33A02C",
          "2006" = "#FB9A99",
          "2008" = "#E31A1C",
          "2011" = "#FDBF6F",
          "2014" = "#FF7F00",
          "2017" =  "#CAB2D6",
          "2021" =  "#6A3D9A"
)



pr1.16s_Prelim2 <- ggplot(data=df.melt2, aes(x=Year, y=value, color = Year)) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA, coef = 1) +  # Remove whiskers) +
  geom_jitter(size=1.5, width = 0.15, alpha = 0.5) +  # Set alpha to 0.5 (50% transparent)) +
  scale_colour_manual(values=cols)+
  theme(axis.text.x  = element_text(angle=90, hjust=1, vjust = 1) ) +
  #facet_wrap("variable",scales = "free", nrow = 1, ncol = 3) +
  labs(x = NULL, y = "Shannon alpha diversity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  guides(alpha = guide_legend(title = "Year Received")) +  # Customize the legend title
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 7.5, by = 2), limits = c(0, 7.5)) +  # Set y-axis breaks
  facet_wrap(~Platform)
pr1.16s_Prelim2

# export data from df
write.csv(df.melt2, file = "/Writing/Prelim study/Extra figures/df_melt2_data.csv", row.names = FALSE)


###Ratios#####
Ratios = read_xlsx("/Writing/Prelim study/Extra figures/Rarified_ratio_LR_divd_SR.xlsx")

#convert year to a factor
Ratios[, 'Year'] <- lapply(Ratios[, 'Year'], factor)
Ratios

library(ggplot2)

#Jitter_plot
ggplot(Ratios, aes(x = Year, y = Ratio_LR_SR)) +
  geom_jitter(width = 0.2, height = 0, size = 3) +
  labs(x = "Year", y = "Ratio LR SR") +
  ggtitle("Jitter Plot of Ratio LR SR by Year")

#Boxplot
ggplot(Ratios, aes(x = Year, y = Ratio_LR_SR)) +
  geom_boxplot() +
  labs(x = "Year", y = "Ratio LR SR") +
  ggtitle("Boxplot of Ratio LR SR by Year")

Rarefied_ratios <- ggplot(Ratios, aes(x = Year, y = Ratio_LR_SR, color = Year)) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers) 
  geom_jitter(aes(color = Year), size = 1.5, width = 0.15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  labs(x = "Year", y = "Shannon alpha diversity") +
  scale_color_manual(values = c(
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Year", y = "Alpha diversity ratio") +
  scale_y_continuous(expand = c(0, 0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), limits = c(0,1.2))  # Set y-axis breaks

Rarefied_ratios


#######Unrarefied#######
phy.Prelim_filtered

number_of_asvs <- ntaxa(phy.Prelim_filtered)

# Print the number of ASVs
print(number_of_asvs)

#####Shannon Diversity Index#########
shan3 <- plot_richness(phy.Prelim_filtered, measures=c("Shannon"))

out3 <- data.frame(sample=shan3$data$id,
                  shannon=shan3$data$value,
                  Year=shan3$data$Year,
                  Platform=shan3$data$Platform)

out3$Eff_No_ASV <- exp(out3$shannon)

names(out3)
# [1] "sample"     "shannon"    "Year"       "Platform"   "Eff_No_ASV"

library(data.table)
df.melt3 <- melt(as.data.table(out3[, c("sample", "Year", "Platform", "shannon")]), id.vars = c("sample", "Year", "Platform"))
df.melt3

str(df.melt3)
# Classes ‘data.table’ and 'data.frame':	172 obs. of  5 variables:
#$ sample  : chr  "AAA5173.fastq" "AAA5228.fastq" "AAA5670.fastq" "AAA5727.fastq" ...
#$ Year    : Factor w/ 10 levels "2001","2002",..: 1 1 1 1 1 1 1 1 3 3 ...
#$ Platform: chr  "PacBio" "PacBio" "PacBio" "PacBio" ...
#$ variable: Factor w/ 1 level "shannon": 1 1 1 1 1 1 1 1 1 1 ...
#$ value   : num  605 413 414 231 508 ...
#- attr(*, ".internal.selfref")=<externalptr>  


unique(df.melt3$Year)
# [1] 2001 2003 2002 2004 2006 2008 2011 2014 2017 2021
#Levels: 2001 2002 2003 2004 2006 2008 2011 2014 2017 2021 

 df.melt3$Year <- factor(df.melt3$Year,
                       levels=c(
                         "2001", "2002", "2003", "2004", "2006", 
                         "2008", "2011", "2014", "2017", "2021"
                       ),
                       labels=c(
                         "2001", "2002", "2003", "2004", "2006", 
                         "2008", "2011", "2014", "2017", "2021"
                       ),
                       ordered = TRUE)

cols <- c("2001" = "#A6CEE3", 
          "2002" = "#1F78B4",
          "2003" = "#B2DF8A",
          "2004" = "#33A02C",
          "2006" = "#FB9A99",
          "2008" = "#E31A1C",
          "2011" = "#FDBF6F",
          "2014" = "#FF7F00",
          "2017" =  "#CAB2D6",
          "2021" =  "#6A3D9A"
)





p16s_Prelim3 <- ggplot(data=df.melt3, aes(x=Year, y=value, color = Year)) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA, coef = 1) +  # Remove whiskers) +
  geom_jitter(size=1.5, width = 0.15, alpha = 0.5) +  # Set alpha to 0.5 (50% transparent)) +
  scale_colour_manual(values=cols)+
  theme(axis.text.x  = element_text(angle=90, hjust=1, vjust = 1) ) +
  #facet_wrap("variable",scales = "free", nrow = 1, ncol = 3) +
  labs(x = NULL, y = "Shannon alpha diversity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  guides(alpha = guide_legend(title = "Year Received")) +  # Customize the legend title
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 7.5, by = 2), limits = c(0, 7.5)) +  # Set y-axis breaks
  facet_wrap(~Platform)
p16s_Prelim3

library(gridExtra)


# Combine plots with labels
combined_plotAlpha_and_ratios <- grid.arrange(
  arrangeGrob(pr1.16s_Prelim2 + theme(legend.position = "none"), top = textGrob("a)", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold"))),
  arrangeGrob(p16s_Prelim3 + theme(legend.position = "none"), top = textGrob("b)", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold"))),
  combined_plotRatios,
  ncol = 1,
  widths = c(1)
)

combined_plotAlpha_and_ratios



p16s_Prelim4 <- ggplot(data=df.melt3, aes(x=Year, y=value, color = Year)) +
  ggtitle("Bacteria alpha diversity Shannon")+
  theme_bw() +
  geom_boxplot() +
  scale_colour_manual(values=cols)+
  theme(axis.text.x  = element_text(angle=90, hjust=1, vjust = 1) ) +
  #facet_wrap("variable",scales = "free", nrow = 1, ncol = 3) +
  labs(x = NULL, y = "Shannon alpha diversity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  guides(alpha = guide_legend(title = "Year Received")) +  # Customize the legend title
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 7.2, by = 2), limits = c(0, 7.5)) +  # Set y-axis breaks
  facet_wrap(~Platform)
p16s_Prelim4


# export data from df
write.csv(df.melt2, file = "/Prelim study/Extra figures/df_melt_unrarified_data.csv", row.names = FALSE)


###Ratios#####
RatiosUR = read_xlsx("/Prelim study/Extra figures/Un_Rarified_ratio_LR_divd_SR.xlsx")

#convert year to a factor
RatiosUR[, 'Year'] <- lapply(RatiosUR[, 'Year'], factor)
RatiosUR


UnRarefied_ratios <- ggplot(RatiosUR, aes(x = Year, y = Ratio, color = Year)) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers)
  geom_jitter(aes(color = Year), size = 1.5, width = 0.15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  scale_color_manual(values = c(
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Year", y = "Alpha diversity ratio") +
  scale_y_continuous(expand = c(0, 0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), limits = c(0,1.2))  # Set y-axis breaks


UnRarefied_ratios

UnRarefied_ratios2 <- ggplot(RatiosUR, aes(x = Year, y = Ratio, color = Year)) +
  theme_bw() +
  geom_boxplot() +  # Add black outlines to the boxes
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  scale_color_manual(values = c(
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Year", y = "Alpha diversity ratio") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2))  # Set y-axis breaks

UnRarefied_ratios2




library(gridExtra)

# Combine plots with labels
combined_plot <- grid.arrange(
  arrangeGrob(Rarefied_ratios + theme(legend.position = "none"), top = textGrob("a)", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold"))),
  arrangeGrob(UnRarefied_ratios + theme(legend.position = "none"), top = textGrob("b)", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold"))),
  ncol = 2,
  widths = c(1, 1)
)

# Print the combined plot
print(combined_plot)

  


###Ratios_rarevsUnRare#####
Ratios_both = read_xlsx("/Writing/Prelim study/Extra figures/RarevsUnRareRatio.xlsx")

#convert year to a factor
Ratios_both[, 'Year'] <- lapply(Ratios_both[, 'Year'], factor)
Ratios_both

pRatios_both <- ggplot(Ratios_both, aes(x = Year, y = Ratio, color = Year)) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA, coef = 1) +  
  geom_jitter(aes(color = Year), size = 1.5, width = 0.15, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  scale_color_manual(values = c(
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666", "#999999", "#E5D8BD")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Year", y = "Alpha diversity ratio") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2)) + # Set y-axis breaks
facet_wrap(~Ratio_Type)

pRatios_both

Ratios_scatter <- ggplot(Ratios_both, aes(x = Year, y = Ratio)) +
  theme_bw() +
  geom_jitter(aes(color = Ratio_Type), size = 1.5, width = 0.15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  scale_color_manual(values = c(
    "#CAB2D6", "#1F78B4")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Year", y = "Alpha diversity ratio") +
  scale_y_continuous(expand = c(0, 0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), limits = c(0,1.2))  # Set y-axis breaks


Ratios_scatter

# Arrange the first two plots in a single column
top_row <- arrangeGrob(
  pr1.16s_Prelim2 + theme(legend.position = "none"), 
  p16s_Prelim3 + theme(legend.position = "none"), 
  ncol = 1
)

# Arrange the second two plots side by side
bottom_row <- arrangeGrob(
  Rarefied_ratios + theme(legend.position = "none"), 
  Ratios_scatter + theme(legend.position = "none"), 
  ncol = 2
)

# Combine the top and bottom rows
combined_plot <- grid.arrange(
  top_row,
  bottom_row,
  ncol = 1,
  heights = c(2, 1)  # Set both rows to have equal height
)

combined_plot

###
# Combine plots with labels
combined_plot <- grid.arrange(
  arrangeGrob(p16s_Prelim3 + theme(legend.position = "none"), top = textGrob("a)", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold"))),
  arrangeGrob(pr1.16s_Prelim2 + theme(legend.position = "none"), top = textGrob("b)", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold"))),
  arrangeGrob(pRatios_both + theme(legend.position = "none"), top = textGrob("c)", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold"))),
  arrangeGrob(Ratios_scatter + theme(legend.position = "right"), top = textGrob("d)", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold"))),
  ncol = 2
)


#########Ratio_scatter_w_shapes_and_colours###
Ratios_scatter <- ggplot(Ratios_both, aes(x = Year, y = Ratio, color = Ratio_Type, shape = Ratio_Type)) +
  theme_bw() +
  geom_jitter(aes(shape = Ratio_Type), size = 3, width = 0.15) +  # Specify shape within geom_jitter
  scale_shape_manual(values = c(16, 15)) +  # 16 for circle, 15 for square
  scale_color_manual(values = c(
    "#CAB2D6", "#1F78B4")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        strip.text = element_text(size = 12, face = 'bold')) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Year", y = "Alpha diversity ratio") +
  scale_y_continuous(expand = c(0, 0), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), limits = c(0,1.2))  # Set y-axis breaks

Ratios_scatter

#########From_PCoA_plot_illumina_PacBio.R#######
library(ggpubr)
# Arrange the plots without the legends
P6 <- ggarrange(UnRarefied_ratios, p1a6, PCoA_Illumina_axis_plot2, PCoA_PacBio_axis_plot2,
                ncol = 2, nrow = 2, 
                labels = c("a)", "b)", "c)", "d)"),
                common.legend = TRUE, 
                widths = 2, heights = 4, 
                align = "hv")

# Print the final plot
P6


########Number of ASVs######
number_of_asvs <- ntaxa(phy.Prelim_filtered)

# Print the number of ASVs
print(number_of_asvs)
#[1] 38572

phy_PacBio <- subset_samples(phy.Prelim_filtered, Platform == "PacBio")
phy_PacBio

min( sample_sums(phy_PacBio) ) # 1563
min(taxa_sums(phy_PacBio)) # 0

# prune taxa that have zero sequence reads
phy_PacBio <- prune_taxa(taxa = taxa_sums(phy_PacBio) > 0, x = phy_PacBio)
phy_PacBio

number_of_asvsPB <- ntaxa(phy_PacBio)

# Print the number of ASVs
print(number_of_asvsPB)
#[1] 14261


# Extract the taxonomic table
tax_table_data <- tax_table(phy_PacBio)

# Extract the genus column (assuming genus is in the 6th column of the taxonomic table)
genera <- tax_table_data[, "Genus"]

# Remove any NA values (optional, depending on your data)
genera <- na.omit(genera)

# Count the unique genera
number_of_genera <- length(unique(genera))

# Print the number of genera
print(number_of_genera)
#[1] 573


# Extract the genus column (assuming genus is in the 6th column of the taxonomic table)
species <- tax_table_data[, "Species"]

# Remove any NA values (optional, depending on your data)
species <- na.omit(species)

# Count the unique species
number_of_species <- length(unique(species))

# Print the number of species
print(number_of_species)
#[1] 1100

otu_table_data <- otu_table(phy_PacBio)
# Calculate the total number of reads for all samples
total_reads <- sum(otu_table_data)
total_reads
#[1] 355749



#####
phy_Illumina <- subset_samples(phy.Prelim_filtered, Platform == "Illumina")
phy_Illumina

min( sample_sums(phy_Illumina) ) # 52271
min(taxa_sums(phy_Illumina)) # 0

# prune taxa that have zero sequence reads
phy_Illumina <- prune_taxa(taxa = taxa_sums(phy_Illumina) > 0, x = phy_Illumina)
phy_Illumina


number_of_asvsIll <- ntaxa(phy_Illumina)

# Print the number of ASVs
print(number_of_asvsIll)
#[1] 24311

# Extract the taxonomic table
tax_table_data2 <- tax_table(phy_Illumina)

# Extract the genus column (assuming genus is in the 6th column of the taxonomic table)
genera2 <- tax_table_data2[, "Genus"]

# Remove any NA values (optional, depending on your data)
genera2 <- na.omit(genera2)

# Count the unique genera
number_of_genera2 <- length(unique(genera2))

# Print the number of genera
print(number_of_genera2)
#[1] 534

# Extract the genus column (assuming genus is in the 6th column of the taxonomic table)
species2 <- tax_table_data2[, "Species"]

# Remove any NA values (optional, depending on your data)
species2 <- na.omit(species2)

# Count the unique species
number_of_species2 <- length(unique(species2))

# Print the number of species
print(number_of_species2)
#[1] 282

# Extract the genus column (assuming genus is in the 6th column of the taxonomic table)
species2 <- tax_table_data2[, "Species"]

# Remove any NA values (optional, depending on your data)
species2 <- na.omit(species2)

# Count the unique species
number_of_species2 <- length(unique(species2))

# Print the number of species
print(number_of_species2)
#[1] 282

otu_table_data2 <- otu_table(phy_Illumina)
# Calculate the total number of reads for all samples
total_reads2 <- sum(otu_table_data2)
total_reads2
#[1] 7340635

# Create a factor corresponding to the Genera
genfac = factor(tax_table(phy_Illumina)[, "Genus"])
# Tabulate the counts for each genera in each sample
gentab = apply(otu_table(phy_Illumina), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
head(gentab)[, 1:10]

  
# To get number of non-zero genera per sample, sum the values that are above
# your threshold, in your case, 1.
observationThreshold = 1
apply(gentab > observationThreshold, 2, sum) 
  
table(tax_table(phy_Illumina)[,"Genus"])
