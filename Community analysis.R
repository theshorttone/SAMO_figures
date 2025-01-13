### Code workthrough source: 
#https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

# initialize --------------------------------------------------------------
library(vegan)
library(phyloseq)
library(ggplot2)
library(dplyr)
setwd("/Users/zacharyshortt/Desktop/Oak_Thesis/bio_informatics/SAMO_FIGS")


# phyloseq import, data organization -------------------------------------------------------------
ps <- readRDS("phyloseq_data_genus_decontam.rds")

# Access the sample_data from the phyloseq object
sample_df <- sample_data(ps)

# Convert to a data.frame if needed
sample_df <- as(sample_df, "data.frame")

# Create the new column 'block_id' based on 'plot_num' and convert it to a factor
sample_df <- sample_df %>%
  mutate(block_id = case_when(
    plot_num %in% c("1a", "1b", "2a", "2b") ~ "Block 1",  
    plot_num %in% c(3, 4) ~ "Block 2",                    
    plot_num %in% c(5, 6) ~ "Block 3",                    
    plot_num %in% c(7, 8) ~ "Block 4",                    
    TRUE ~ NA_character_                                  
  ))

# Convert 'block_id' to a factor explicitly
sample_df$block_id <- factor(sample_df$block_id)

# Reassign the modified sample_data back to the phyloseq object
sample_data(ps) <- sample_df

####get rid of controls
ps_new <- subset_samples(ps, sample_id != "neg-control1" & sample_id != "neg-control2")


# combine samples  --------------------------------------------------------sample_data(ps)$Tree.ID
ps_combined <- merge_samples(ps_new, "Tree.ID")

# Summarize metadata for each Tree.ID
metadata <- sample_data(ps_new) %>% as.data.frame()
new_metadata <- metadata %>%
  group_by(Tree.ID) %>%
  summarise(across(everything(), ~ if (is.numeric(.)) mean(.) else first(.))) %>%
  as.data.frame() # Convert tibble to data frame

# Ensure row names match sample names in ps_combined
row.names(new_metadata) <- new_metadata$Tree.ID

# Assign back the updated sample data
sample_names(ps_combined) <- new_metadata$Tree.ID  # Ensure names align
sample_data(ps_combined) <- sample_data(new_metadata)

# verification that combo worked -----------------------------------------------------------------------

# Extract Tree.ID from original ps
original_tree_ids <- sample_data(ps_new)$Tree.ID

# Extract Tree.ID from combined ps_combined
combined_tree_ids <- row.names(sample_data(ps_combined))

all(combined_tree_ids %in% original_tree_ids) # Should return TRUE

# Check original metadata
original_metadata <- sample_data(ps_new) %>% as.data.frame()

# Check combined metadata
combined_metadata <- sample_data(ps_combined) %>% as.data.frame()

# Join and compare
comparison <- dplyr::full_join(
  original_metadata,
  combined_metadata,
  by = "Tree.ID"
)

# View the comparison
View(comparison) # Or print(comparison) for console

otu_matrix_og <- otu_table(ps_new)

otu_matrix_com <- otu_table(ps_combined)

#  family, phylum dataframe creation----------------------------------------------------------------------

##############READ
##in the future it may be best to do the ps family and ps phylum stuff before you rarefy to take out samples, and just rarefy for the ndms and alpha diversity below 

ps_phylum <- ps_combined %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01)                       # Filter out low abundance taxa

ps_family_05 <- ps_combined %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.005) 

ps_family_20 <- ps_combined %>%
  tax_glom(taxrank = "Family") %>%                     # Agglomerate at Family level
  transform_sample_counts(function(x) {x / sum(x)}) %>% # Transform to relative abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.005) %>%                        # Filter out very low abundance values
  mutate(Family = ifelse(Abundance < 0.2, "Other", Family)) # Group families < 20% into "Other"


# Set colors for plotting
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)


# plotting family richness by innoc------------------------------------------------


ggplot(ps_family_05 %>% filter(tree.species != "agrifolia " & tree.species != "berbifolia/agrifolia" & tree.species != "pine"), 
       aes(x = tree.species, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +                    # Corrected by adding "+" to chain layers
  theme(legend.position = "right",                 # Position the legend to the right
        legend.title = element_text(size = 8),      # Smaller legend title
        legend.text = element_text(size = 6),       # Smaller legend text
        legend.key.size = unit(0.5, "lines"),       # Smaller legend key size (color box)
        legend.spacing = unit(0.5, "lines")) +      # Adjust space between legend items
  facet_wrap(~colonized)                           # Facet by the 'colonized' variable

ggsave("Figures/all_familyrichness_bycol.png")

ggplot(ps_family_20 %>% filter(tree.species != "agrifolia " & tree.species != "berbifolia/agrifolia" & tree.species != "pine"), 
       aes(x = tree.species, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +                    # Corrected by adding "+" to chain layers
  theme(legend.position = "right",                 # Position the legend to the right
        legend.title = element_text(size = 8),      # Smaller legend title
        legend.text = element_text(size = 6),       # Smaller legend text
        legend.key.size = unit(0.5, "lines"),       # Smaller legend key size (color box)
        legend.spacing = unit(0.5, "lines")) +      # Adjust space between legend items
  facet_wrap(~colonized)                           # Facet by the 'colonized' variable

ggsave("Figures/common_familyrichness_bycol.png")

# plotting family richness by species------------------------------------------------

ggplot(ps_family_05 %>% filter(tree.species != "agrifolia " & tree.species != "berbifolia/agrifolia"), 
       aes(x = tree.species, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  
  theme(legend.position = "right",                 # Position the legend to the right
      legend.title = element_text(size = 8),      # Smaller legend title
      legend.text = element_text(size = 6),       # Smaller legend text
      legend.key.size = unit(0.5, "lines"),       # Smaller legend key size (color box)
      legend.spacing = unit(0.5, "lines"))

ggsave("Figures/all_familyrichness_pine.png")


ggplot(ps_family_20 %>% filter(tree.species != "agrifolia " & tree.species != "berbifolia/agrifolia"), 
       aes(x = tree.species, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") 
  #scale_fill_manual(values = phylum_colors)
ggsave("Figures/common_familyrichness.png")
 
# plotting family richness by plot------------------------------------------------

ggplot(ps_family_20, aes(x = plot_num, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") 
  #scale_fill_manual(values = phylum_colors)
ggsave("Figures/common_familyrichness_byplot.png")


# plotting phylum richness ------------------------------------------------
ggplot(ps_phylum %>% filter(tree.species != "agrifolia " & tree.species != "berbifolia/agrifolia"), 
       aes(x = tree.species, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors)

ggsave("Figures/common_phylumrichness_byspecie.png")


ggplot(ps_phylum, aes(x = block_id, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors)

ggsave("Figures/common_phylumrichness_byblock.png")

ggplot(ps_phylum %>% filter(tree.species != "agrifolia " & tree.species != "berbifolia/agrifolia" & tree.species != "pine"), 
       aes(x = plot_num, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  facet_wrap(~tree.species)  # Facet by species
ggsave("Figures/common_phylumrichness_by_plot_species.png")


# rarefaction curve -------------------------------------------------------

set.seed(711)
rarefied_ps <-rarefy_even_depth(ps_combined, sample.size = min(sample_sums(ps_combined)),
                                rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
# Histogram of sample read counts

# Make a data frame with a column for the read counts of each sample

otu_matrix <- as(otu_table(ps_combined), "matrix")
otu_matrix_t <- t(otu_matrix)
rarecurve(otu_matrix, step=50, cex=0.7)

ggsave("Figures/rarifaction_curve.png")


sample_sizes <- rowSums(otu_matrix)

# Identify the sample(s) with the smallest size
sample_sizes <- rowSums(otu_matrix)

# Identify the 5 samples with the smallest total sizes
sorted_samples <- sort(sample_sizes) # Sort samples by size
smallest_5 <- head(sorted_samples, 5) # Get the smallest 5
smallest_5_samples <- names(smallest_5)
cat("Samples with the smallest 5 total sizes:\n")
print(smallest_5)


sample_sum_df <- data.frame(sum = sample_sums(ps_combined))

ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())


# NDMS --------------------------------------------------------------------


# Ordinate generally to compare more variables than what I'm using 
rarefied_ndms<- ordinate(
  physeq = rarefied_ps, 
  method = "NMDS", 
  distance = "bray"
)
# Plot 
plot_ordination(
  physeq = rarefied_ps,
  ordination = rarefied_ndms,
  color = "plot_num",
  shape = "tree.species",
  title = "NMDS of Fungal communities"
)

ggsave("Figures/NDMS.png")

# NDMS stats------------------------------------------------------


set.seed(1)

# Calculate bray curtis distance matrix
bray <- phyloseq::distance(rarefied_ps, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(rarefied_ps))

# Adonis test
adonis(bray ~ plot_num, data = sampledf)

beta <- betadisper(bray, sampledf$plot_num)
permutest(beta)

###########################################################################################


####Constrained Ordinates 
ps_not_na <- rarefied_ps %>%
  subset_samples(
    !is.na(tree.species) & 
      tree.species != "agrifolia " & 
      tree.species != "berbifolia/agrifolia" & # Exclude entries where tree.species is "agrifolia"
      !is.na(block_id)  # Keep other filtering conditions
  )

ps_not_na_dis <- phyloseq::distance(physeq = ps_not_na, method = "bray")


# CAP ordinate
cap_ord <- ordinate(
  physeq = ps_not_na, 
  method = "CAP",
  distance = ps_not_na_dis,
  formula = ~ tree.species + block_id
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = ps_not_na, 
  ordination = cap_ord, 
  color = "block_id", 
  axes = c(1,2)
) + 
  aes(shape = tree.species) + 
  geom_point(aes(colour = block_id), alpha = 0.8, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) 
 # scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", 
                         #       "#1919ff", "darkorchid3", "magenta", "orange", "black", "light"))
  


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    data = arrowdf, 
    show.legend = FALSE
  )
ggsave("Figures/CAP.png")

###########################################################################################

###############################################################################

###########################################################################################
###ALPHA DIVERSITY. NOT THAT HELPFUL, the family level stuff was better I think 
# Initialize matrices to store richness and evenness estimates
min_lib <- min(sample_sums(ps_combined))


nsamp = nsamples(ps_combined)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(ps_combined)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(ps_combined)

# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(ps_combined, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
  }

# Create a new dataframe to hold the means and standard deviations of richness estimates
sample_id <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(sample_id, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
sample_id <- row.names(evenness)
mean_even <- apply(evenness, 1, mean)
sd_even <- apply(evenness, 1, sd)
measure_even <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(sample_id, mean_even, sd_even, measure_even)

alpha <- merge(rich_stats, even_stats, by = "sample_id")

s <- data.frame(sample_data(ps_combined))

# Convert row names to a column in s
s$sample_id <- rownames(s)

alphadiv <- merge(alpha, s, by = "sample_id") 


ggplot(alphadiv, aes(x = block_id, y = mean, color = tree.species, group = tree.species, shape = tree.species)) +
  geom_point(size = 2) + 
  geom_line(size = 0.8) +
  #facet_wrap(~measure, ncol = 1, scales = "free") +
  #scale_color_manual(values = c("#E96446", "#302F3D", "#87CEFA")) +
  scale_x_discrete(
    drop = FALSE
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

sum_alpha <- alphadiv%>%
  group_by(tree.species, plot_num, colonized) %>%
    summarise(
      mean_div = mean(mean_even, na.rm = TRUE),
      sd_div = sd(mean_even, na.rm = TRUE),
      mean_rich = mean(mean, na.rm = TRUE),
      sd_rich = sd(mean, na.rm = TRUE),
      n = n()
    ) %>%
  mutate(
    se_div = sd_div / sqrt(n),
    se_rich = sd_rich / sqrt(n)
  )

ggplot(sum_alpha, aes(x = plot_num, y = mean_rich, fill = tree.species, group = tree.species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +  # Adjust dodge and bar width
  geom_errorbar(
    aes(ymin = mean_rich- se_rich, ymax = mean_rich + se_rich),
    position = position_dodge(width = 0.8),  # Ensure error bars align with bars
    width = 0.25  # Width of the error bars
  ) 
ggsave("Figures/mean_rich.png")

ggplot(sum_alpha, aes(x = plot_num, y = mean_div, fill = tree.species, group = tree.species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +  # Adjust dodge and bar width
  geom_errorbar(
    aes(ymin = mean_div- se_div, ymax = mean_div + se_div),
    position = position_dodge(width = 0.8),  # Ensure error bars align with bars
    width = 0.25  # Width of the error bars
  ) 

ggsave("Figures/mean_div.png")

