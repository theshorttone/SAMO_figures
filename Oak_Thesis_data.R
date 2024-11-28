library(readxl)
library(tidyverse)
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
# purple_pal <- colorRampPalette(brewer.pal("Purples"))
# orange_pal <-colorRampPalette(brewer.pal("Oranges"))
colonization_data <- read_excel("/Users/zacharyshortt/Desktop/DNA_work_SMM/Colonization_data.xlsx")
master_data <-read_excel("/Volumes/f004qd8/Samo_Data/master_data.xlsx")

combined_coldata <- colonization_data %>%
  group_by(`Tree code`) %>%
  summarize(across(everything(), sum, na.rm = TRUE), .groups = 'drop')
#print(combined_coldata)

# Merging the data frames
merged_data_all <- merge(master_data,combined_coldata , by = "Tree code", all = TRUE)
#print(combined_coldata)
#colnames(merged_data_all)

# rename columns with annoying names
merged_data_all <- merged_data_all %>% rename(survival = `Alive (Y/N)`, col = `Total col`, uncol = `uncolonized tips (with E)`)

#calculating %colonization
merged_data_all <- merged_data_all %>%
  mutate(
   percent_col = (col/(uncol+col))*100 ##using EG data since E is definatly root hairs, and this excludes E, might want to revist and add E and G as uncolonized
  )

##Changing height data to cm
merged_data_all <- merged_data_all %>%
  mutate(
    Heightin_6_8_24 = Heightin_6_8_24 * 2.54,
    Heightin_12_18_23 = Heightin_12_18_23 * 2.54# Convert height to cm
  )

# Calculate radius and volume (cone model)
merged_data_all <- merged_data_all %>%
  mutate(
    radius = width_7_5_24 / 2,  # Radius in cm (width is already in cm)
    volume = (1/3) * pi * (radius^2) * Heightin_6_8_24  # Volume calculation
  )

# Plot the volume against species using ggplot2
ggplot(merged_data_all, aes(x = species, y = volume, fill = species)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Tree Volume by Species",
    x = "Species",
    y = "Volume (cm³)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Plot volume by species and innoculation status 
# Step 1: Calculate mean and standard error for each species and inoculation status
summary_data <- merged_data_all %>%
  group_by(species, innoculation) %>%
  summarise(
    mean_volume = mean(volume, na.rm = TRUE),  # Mean volume
    sd_volume = sd(volume, na.rm = TRUE),      # Standard deviation
    n = n()                                    # Number of observations
  ) %>%
  mutate(
    se_volume = sd_volume / sqrt(n)  # Calculate standard error (SE)
  )

# Plot with side-by-side bars and error bars
ggplot(summary_data, aes(x = species, y = mean_volume, fill = innoculation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +  # Adjust dodge and bar width
  geom_errorbar(
    aes(ymin = mean_volume - se_volume, ymax = mean_volume + se_volume),
    position = position_dodge(width = 0.8),  # Ensure error bars align with bars
    width = 0.25  # Width of the error bars
  ) +
  labs(
    title = "Tree Volume by Species and Inoculation Status with Error Bars",
    x = "Species",
    y = "Mean Volume (cm³)",
    fill = "Inoculation Status"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


summary_data_col <- merged_data_all %>%
  filter(!is.na(percent_col)) %>%            # Remove rows with NA in percent_col
  group_by(species, innoculation) %>%
  summarise(
    col_avg = mean(percent_col, na.rm = TRUE),  # Mean colonization
    sd_col = sd(percent_col, na.rm = TRUE),     # Standard deviation
    n = n()                                     # Number of observations
  ) %>%
  mutate(
    se_col = sd_col / sqrt(n - 1)  # Calculate standard error (SE)
  )


ggplot(summary_data_col, aes(x = species, y = col_avg, fill = innoculation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +  # Adjust dodge and bar width
  geom_errorbar(
    aes(ymin = col_avg - se_col, ymax = col_avg + se_col),
    position = position_dodge(width = 0.8),  # Ensure error bars align with bars
    width = 0.25  # Width of the error bars
  ) +
  labs(
    title = "Tree colonization by Species and Inoculation Status with Error Bars",
    x = "Species",
    y = "Mean Colonization (%)",
    fill = "Inoculation Status"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


#looking at plots 
summary_data_col_plot <- merged_data_all %>%
  filter(!is.na(percent_col)) %>%            # Remove rows with NA in percent_col
  group_by(species, innoculation, plot_num) %>%
  summarise(
    col_avg = mean(percent_col, na.rm = TRUE),  # Mean colonization
    sd_col = sd(percent_col, na.rm = TRUE),     # Standard deviation
    n = n()                                     # Number of observations
  ) %>%
  mutate(
    se_col = sd_col / sqrt(n - 1)  # Calculate standard error (SE)
  )


# ggplot(summary_data_col_plot, aes(x = species, y = col_avg, fill = plot_num)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
#   geom_errorbar(
#     aes(ymin = col_avg - se_col, ymax = col_avg + se_col),
#     position = position_dodge(width = 0.8),
#     width = 0.25
#   ) +
#   labs(
#     title = "Tree colonization by Species and Inoculation Status with Error Bars",
#     x = "Species",
#     y = "Mean Colonization (%)",
#     fill = "Plot Number"  # Adjust legend title
#   ) +
#   scale_fill_manual(values = c("1b" = "grey", "2b" = "grey", "2a" = "grey", "3"= "blue", "4" = "red", "5" = "red", "6" = "blue", "7"= "red", "8" = "blue")) +  # Customize colors
#   theme_minimal() +
#   theme(legend.position = "top")

# Filter out plot numbers 1b, 2a, and 2b
summary_data_col_plot_filtered <- summary_data_col_plot %>%
  filter(!plot_num %in% c("1b", "2a", "2b"))



# Define colors
colors_purples <- brewer.pal(9, "Purples")[5:7]  # Medium-dark shades of Purples
colors_oranges <- brewer.pal(9, "Oranges")[5:7]  # Medium-dark shades of Oranges
custom_colors <- c(colors_purples, colors_oranges)

# Define mapping of plot_num to colors
color_mapping <- c("3" = custom_colors[1], 
                   "6" = custom_colors[2], 
                   "8" = custom_colors[3],
                   "4" = custom_colors[4],
                   "5" = custom_colors[5],
                   "7" = custom_colors[6])

# Ensure plot_num is a factor with correct order
summary_data_col_plot_filtered <- summary_data_col_plot_filtered %>%
  mutate(plot_num = factor(plot_num, levels = c(3, 6, 8, 4, 5, 7)))

# Plot
ggplot(summary_data_col_plot_filtered, aes(x = species, y = col_avg, fill = plot_num)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_errorbar(
    aes(ymin = col_avg - se_col, ymax = col_avg + se_col),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  geom_text(
    aes(label = plot_num, y = col_avg + se_col + 2),  # Position text above the error bar
    position = position_dodge(width = 0.8),
    vjust = 0,
    size = 4,
    color = "black"
  ) +
  labs(
    title = "Tree Colonization by Species and Inoculation Status with Error Bars",
    x = "Species",
    y = "Mean Colonization (%)",
    fill = "Plot"
  ) +
  scale_fill_manual(values = color_mapping)

##############


ggplot(summary_data_col, aes(x = species, y = col_avg, fill = innoculation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +  # Adjust dodge and bar width
  geom_errorbar(
    aes(ymin = col_avg - se_col, ymax = col_avg + se_col),
    position = position_dodge(width = 0.8),  # Ensure error bars align with bars
    width = 0.25  # Width of the error bars
  ) +
  labs(
    #title = "",
    x = "Species",
    y = "Mean Volume (cm³)",
    fill = "Inoculation Status"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

##Now all the data should be downloaded, below is graph creation 


ggplot(merged_data_all, aes(x = species, fill = survival)) +
  geom_bar(position = "dodge") +
  labs(title = "Survival by Tree Species",
       x = "Tree Species",
       y = "Count") +
  theme_minimal()
# Sample data frame

##calculating percent survival for each species, and corresponding chi squared 

# Create a contingency table
contingency_table <- table(merged_data_all$species, merged_data_all$survival)

# Perform chi-squared test
chi_squared_test <- chisq.test(contingency_table)

print(chi_squared_test)

###

# Calculate percent survival and standard error
percent_survival <- merged_data_all %>%
  group_by(species) %>%
  summarise(
    survival_rate = mean(survival == "Y") * 100,
    n = n(),  # Number of observations
    se = sqrt((survival_rate / 100) * (1 - survival_rate / 100) / n) * 100  # Standard Error
  )

# View the calculated percent survival and standard error
print(percent_survival)

# Create a bar plot for percent survival with error bars
ggplot(percent_survival, aes(x = species, y = survival_rate, fill = species)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = survival_rate - se, ymax = survival_rate + se), width = 0.2) +
  labs(title = "Percent Survival by Species",
       x = "Species",
       y = "Percent Survival (%)") +
  theme_minimal()



## surival by colonizated or not:

contingency_table <- table(merged_data_all$species, merged_data_all$survival)

# Perform chi-squared test
chi_squared_test <- chisq.test(contingency_table)

print(chi_squared_test)

# Calculate percent survival and standard error
percent_survival <- merged_data_all %>%
  group_by(innoculation) %>%
  summarise(
    survival_rate = mean(survival == "Y") * 100,
    n = n(),  # Number of observations
    se = sqrt((survival_rate / 100) * (1 - survival_rate / 100) / n) * 100  # Standard Error
  )

# View the calculated percent survival and standard error
print(percent_survival)

# Create a bar plot for percent survival with error bars
ggplot(percent_survival, aes(x = innoculation, y = survival_rate, fill = innoculation)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = survival_rate - se, ymax = survival_rate + se), width = 0.2) +
  labs(title = "Percent Survival by Species",
       x = "Species",
       y = "Percent Survival (%)") +
  theme_minimal()



library(lme4)
library(lmerTest)

merged_data_all <- merged_data_all %>%
  mutate(block_id = case_when(
    plot_num %in% c("1a", "1b", "2a", "2b") ~ "Block 1",  # If plot_num is 1 or 2, assign "Block 1"
    plot_num %in% c(3, 4) ~ "Block 2",  # If plot_num is 3 or 4, assign "Block 2"
    plot_num %in% c(5, 6) ~ "Block 3",  # If plot_num is 1 or 2, assign "Block 1"
    plot_num %in% c(7, 8) ~ "Block 4",  # If plot_num is 3 or 4, assign "Block 2"
    TRUE ~ NA_character_               # In case there are other values of plot_num, assign NA
  ))


model = lmer(percent_col ~ innoculation * species + (1 | plot_num), data = merged_data_all)
summary(model)

model_block = lmer(percent_col ~ innoculation * species + (1 | block_id/plot_num), data = merged_data_all)

summary(model_block)

model_block2 = lmer(percent_col ~ innoculation * species + (1 | block_id), data = merged_data_all)

summary(model_block2)

model = lmer(volume ~ innoculation * species + (1 | plot_num), data = merged_data_all)
summary(model)

model_block = lmer(volume ~ innoculation * species + (1 | block_id/plot_num), data = merged_data_all)

summary(model_block)

model_block2 = lmer(volume ~ innoculation * species + (1 | block_id), data = merged_data_all)

summary(model_block2)
