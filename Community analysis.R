### Code workthrough source: 
#https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

library(vegan)
library(phyloseq)
library(ggplot2)
library(dplyr)

setwd("/Users/zacharyshortt/Desktop/Oak_Thesis/SAMO_FIGS")

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

library(phyloseq)
library(dplyr)

# Assuming 'physeq' is your phyloseq object

# Step 1: Extract sample data and OTU table from phyloseq
sample_data_df <- as.data.frame(sample_data(ps_new))
otu_table_df <- as.data.frame(otu_table(ps_new))

# Step 2: Aggregate by 'tree_code'
# For categorical columns, you can use 'first' or 'last' to keep a single value per group
# For numerical columns, we will sum them

aggregated_sample_data <- sample_data_df %>%
  group_by(Tree.ID) %>%
  summarise(across(everything(), ~ first(.)))  # Keep the first categorical value for each group

# Sum the OTU table by 'tree_code'
aggregated_otu_table <- otu_table_df %>%
  group_by(tree_code) %>%
  summarise(across(everything(), ~ sum(.)))  # Sum the OTU counts

# Step 3: Recreate the phyloseq object with the aggregated data
# Convert back to appropriate phyloseq components
aggregated_sample_data_ps <- sample_data(aggregated_sample_data)
aggregated_otu_table_ps <- otu_table(as.matrix(aggregated_otu_table), taxa_are_rows = TRUE)

# Recreate the phyloseq object with the aggregated data
aggregated_physeq <- phyloseq(
  aggregated_otu_table_ps,
  aggregated_sample_data_ps,
  tax_table(physeq),  # Retain original tax_table (no change needed here)
  phy_tree(physeq)    # Retain original phylogenetic tree (if applicable)
)

# Check the result
aggregated_physeq




# system.time(rarecurve(otu_table(ps), step = 300, col = "blue", label = FALSE))
##?set.seed

##############READ
##in the future it may be best to do the ps family and ps phylum stuff before you rarefy to take out samples, and just rarefy for the ndms and alpha diversity below 

set.seed(711)
rarefied_ps <-rarefy_even_depth(ps_new, sample.size = min(sample_sums(ps_new)),
                  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# ###RARECURVE
# otu_matrix <- as(otu_table(rarefied_ps), "matrix")
# otu_matrix_t <- t(otu_matrix)
# rarecurve(otu_matrix, step=50, cex=0.1)


# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(ps))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

ps_phylum <- ps_new %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01)                       # Filter out low abundance taxa

ps_family <- ps_new %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.005) 

ps_family_10 <- ps_new %>%
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


ggplot(ps_family %>% filter(tree.species != "agrifolia " & tree.species != "berbifolia/agrifolia" & tree.species != "pine"), 
       aes(x = tree.species, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +                    # Corrected by adding "+" to chain layers
  theme(legend.position = "right",                 # Position the legend to the right
        legend.title = element_text(size = 8),      # Smaller legend title
        legend.text = element_text(size = 6),       # Smaller legend text
        legend.key.size = unit(0.5, "lines"),       # Smaller legend key size (color box)
        legend.spacing = unit(0.5, "lines")) +      # Adjust space between legend items
  facet_wrap(~colonized)                           # Facet by the 'colonized' variable

      

ggplot(ps_family %>% filter(tree.species != "agrifolia" & tree.species != "berbifolia/agrifolia"), 
       aes(x = tree.species, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +


ggplot(ps_family_10 %>% filter(tree.species != "agrifolia " & tree.species != "berbifolia/agrifolia"), aes(x = tree.species, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") 
  #scale_fill_manual(values = phylum_colors)
ggplot(ps_family %>% filter(tree.species != "agrifolia" & tree.species != "berbifolia/agrifolia"), 
       aes(x = tree.species, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
 

ggplot(ps_family_10, aes(x = plot_num, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") 
  #scale_fill_manual(values = phylum_colors)

##PHYLUM####
ggplot(ps_phylum %>% filter(tree.species != "agrifolia " & tree.species != "berbifolia/agrifolia"), aes(x = tree.species, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors)
 

ggplot(ps_phylum, aes(x = block_id, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors)


ggplot(ps_phylum %>% filter(tree.species != "agrifolia " & tree.species != "berbifolia/agrifolia" & tree.species != "pine"), aes(x = plot_num, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  facet_wrap(~tree.species)  # Facet by species

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
set.seed(1)

# Calculate bray curtis distance matrix
bray <- phyloseq::distance(rarefied_ps, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(rarefied_ps))

# Adonis test
adonis(bray ~ plot_num, data = sampledf)

beta <- betadisper(bray, sampledf$plot_num)
permutest(beta)

###############################################################################

###########################################################################################

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
###########################################################################################

###############################################################################

###########################################################################################
###ALPHA DIVERSITY. NOT THAT HELPFUL, the family level stuff was better I think 
# Initialize matrices to store richness and evenness estimates
min_lib <- min(sample_sums(ps_new))


nsamp = nsamples(ps_new)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(ps_new)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(ps_new)

# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(ps_new, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
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

s <- data.frame(sample_data(ps_new))

# Convert row names to a column in s
s$sample_id <- rownames(s)

alphadiv <- merge(alpha, s, by = "sample_id") 


ggplot(alphadiv, aes(x = block_id, y = mean, color = tree.species, group = tree.species, shape = tree.species)) +
  geom_point(size = 2) + 
  geom_line(size = 0.8) +
  facet_wrap(~measure, ncol = 1, scales = "free") +
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

ggplot(sum_alpha, aes(x = plot_num, y = mean_div, fill = tree.species, group = tree.species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +  # Adjust dodge and bar width
  geom_errorbar(
    aes(ymin = mean_div- se_div, ymax = mean_div + se_div),
    position = position_dodge(width = 0.8),  # Ensure error bars align with bars
    width = 0.25  # Width of the error bars
  ) 


