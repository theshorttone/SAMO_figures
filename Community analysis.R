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
# system.time(rarecurve(otu_table(ps), step = 300, col = "blue", label = FALSE))
##?set.seed
set.seed(711)
rarefied_ps <-rarefy_even_depth(ps_new, sample.size = min(sample_sums(ps_new)),
                  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# ###RARECURVE
# otu_matrix <- as(otu_table(rarefied_ps), "matrix")
# otu_matrix_t <- t(otu_matrix)
# rarecurve(otu_matrix, step=50, cex=0.1)

ps<-rarefied_ps
# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(ps))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

ps_phylum <- ps %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01)                       # Filter out low abundance taxa

ps_family <- ps %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.005) 

ps_family_10 <- ps %>%
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





# Ordinate
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
  title = "NMDS of Lake Erie bacterial Communities"
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


