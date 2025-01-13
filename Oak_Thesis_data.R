
# intro -------------------------------------------------------------------

library(readxl)
library(tidyverse)
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

# data call ---------------------------------------------------------------

colonization_data <- read_excel("/Users/zacharyshortt/Desktop/DNA_work_SMM/Colonization_data.xlsx")
master_data <-read_excel("/Users/zacharyshortt/Desktop/DNA_work_SMM/Master_dataaa.xlsx")



# data manipulation, column renaming, ect ---------------------------------

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
    Heightin_12_18_23 = Heightin_12_18_23 * 2.54,# Convert height to cm
    Heightin_12_2024 = Heightin_12_18_23 * 2.54,
    width_7_5_24 = width_7_5_24/10, #convert width to cm,
    width_12_2024 = width_12_2024/10 #convert width to cm
    
  )

# Calculate radius and volume (cone model)
merged_data_all <- merged_data_all %>%
  mutate(
    radius = width_12_2024 / 2,  # Radius in mm (width is already in mm)
    volume = (1/3) * pi * (radius^2) * Heightin_12_2024  # Volume calculation
  )



# Volume bar vs. species graph----------------------------------------------------------------

ggplot(merged_data_all, aes(x = species, y = volume, fill = species)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Tree Volume by Species",
    x = "Species",
    y = "Volume (cm³)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("Figures.mean_species_volume")

# volume, species by innoc bar chart with error bars ---------------------------------------------------------

ggplot(summary_data_growth, aes(x = species, y = mean_volume, fill = innoculation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +  # Adjust dodge and bar width
  geom_errorbar(  aes(ymin = mean_volume - se_volume, ymax = mean_volume + se_volume),
                  position = position_dodge(width = 0.8), width = 0.6)+
  labs(
    x = "Species",
    y = "Mean Volume (cm³)"
  )

# Calculating mean and standard error for each species and inoculation status
summary_data_growth <- merged_data_all %>%
  filter(!is.na(Heightin_12_2024)) %>%  # Ensure no NAs in height variable
  group_by(species, innoculation) %>%
  summarise(
    mean_volume = mean(volume, na.rm = TRUE),  # Mean volume
    sd_volume = sd(volume, na.rm = TRUE),      # Standard deviation of volume
    mean_height = mean(Heightin_12_2024, na.rm = TRUE),  # Mean height
    sd_height = sd(Heightin_12_2024, na.rm = TRUE),      # Standard deviation of height
    n = n()  # Number of observations
  ) %>%
  mutate(
    se_volume = sd_volume / sqrt(n),  # Calculate standard error for volume
    se_height = sd_height / sqrt(n)   # Calculate standard error for height
  )

# Plot volume, species, color by innoc
ggplot(summary_data_growth, aes(x = species, y = mean_volume, fill = innoculation)) +
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



# hight, species by innoc bar chart with error bars-------------------------------------------------

# Plot height, species, color by innoc
ggplot(summary_data_growth, aes(x = species, y = mean_height, fill = innoculation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +  # Adjust dodge and bar width
  geom_errorbar(
    aes(ymin = mean_height - se_height, ymax = mean_height + se_height),
    position = position_dodge(width = 0.8),  # Ensure error bars align with bars
    width = 0.25  # Width of the error bars
  ) +
  labs(
    x = "Species",
    y = "Mean Height (cm)",
    fill = "Inoculation Status"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("Figures.species_mean_ height")

# sumamry data for height by plot------------------------------------------------

summary_data_growth_plot <- merged_data_all %>%
  filter(!is.na(Heightin_12_2024)) %>%            # Remove rows with no height measurements 
  group_by(species, innoculation, plot_num) %>%
  summarise(
    mean_volume = mean(volume, na.rm = TRUE),  # Mean volume
    sd_volume = sd(volume, na.rm = TRUE),      # Standard deviation of volume
    mean_height = mean(Heightin_12_2024, na.rm = TRUE),  # Mean height
    sd_height = sd(Heightin_12_2024, na.rm = TRUE),      # Standard deviation of height
    n = n()  # Number of observations
  ) %>%
  mutate(
    se_volume = sd_volume / sqrt(n),  # Calculate standard error for volume
    se_height = sd_height / sqrt(n)   # Calculate standard error for height
  )


# color scheme for figures  -----------------------------------------------
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
summary_data_growth_plot <- summary_data_growth_plot  %>%
  mutate(plot_num = factor(plot_num, levels = c(3, 6, 8, 4, 5, 7)))



# plotting height, plot species  ------------------------------------------

ggplot(summary_data_growth_plot, aes(x = species, y = mean_height, fill = plot_num)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_height - se_height, ymax = mean_height + se_height),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  geom_text(
    aes(label = plot_num, y = mean_height + se_height + 2),  # Position text above the error bar
    position = position_dodge(width = 0.8),
    vjust = 0,
    size = 4,
    color = "black"
  ) +
  labs(
    title = "Tree Height by Species and Inoculation Status with Error Bars",
    x = "Species",
    y = "Mean Colonization (%)",
    fill = "Plot"
  ) +
  scale_fill_manual(values = color_mapping, labels = 
                      c("3" = "3 innoculated" ,
                        "6" = "6 innoculated" ,
                        "8" = "8 innoculated" ,
                        "4" = "4 control" ,
                        "5"= "5 control" ,
                        "7"= "7 control")
  )

# Plot volumne 
ggplot(summary_data_growth_plot, aes(x = species, y = mean_volume, fill = plot_num)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_volume - se_volume, ymax = mean_volume + se_volume),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  geom_text(
    aes(label = plot_num, y = mean_volume + se_volume + 2),  # Position text above the error bar
    position = position_dodge(width = 0.8),
    vjust = 0,
    size = 4,
    color = "black"
  ) +
  labs(
    title = "Tree Volume by Species and Inoculation Status with Error Bars",
    x = "Species",
    y = "Mean Volume (cm³)",
    fill = "Plot"
  ) +
  scale_fill_manual(values = color_mapping, labels = 
                      c("3" = "3 innoculated" ,
                        "6" = "6 innoculated" ,
                        "8" = "8 innoculated" ,
                        "4" = "4 control" ,
                        "5"= "5 control" ,
                        "7"= "7 control")
  )
###


# removing NA rows for colonization  -------------------------------------------


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


# colonization across innoc and species  ---------------------------------------------------------------


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
    y = "Mean Height(cm)",
    fill = "Inoculation Status"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

summary_data_col_plot <- merged_data_all %>%
  filter(!is.na(percent_col)) %>%            # Remove rows with NA in percent_col
  group_by(species, innoculation, plot_num) %>%
  summarise(

# more data management, remove the top plots ---------------------------------------------------

    
    col_avg = mean(percent_col, na.rm = TRUE),  # Mean colonization
    sd_col = sd(percent_col, na.rm = TRUE),     # Standard deviation
    n = n()                                     # Number of observations
  ) %>%
  mutate(
    se_col = sd_col / sqrt(n - 1)  # Calculate standard error (SE)
  )

summary_data_col_plot_filtered <- summary_data_col_plot %>%
  filter(!plot_num %in% c("1b", "2a", "2b"))

# Ensure plot_num is a factor with correct order
summary_data_col_plot_filtered <- summary_data_col_plot_filtered %>%
  mutate(plot_num = factor(plot_num, levels = c(3, 6, 8, 4, 5, 7)))

# plotting colonization by plots  -----------------------------------------


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
  scale_fill_manual(values = color_mapping, labels = 
                      c("3" = "3 innoculated" ,
                       "6" = "6 innoculated" ,
                       "8" = "8 innoculated" ,
                       "4" = "4 control" ,
                       "5"= "5 control" ,
                       "7"= "7 control")
                    )


# broken idk  -----------------------------------------


ggplot(summary_data, aes(x = species, y = mean_volume, fill = innoculation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +  # Adjust dodge and bar width
  geom_errorbar(
    aes(ymin = mean_volume - se_volume, ymax = mean_volume + se_volume),
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



ggplot(merged_data_all, aes(x = species, fill = survival)) +
  geom_bar(position = "dodge") +
  labs(title = "Survival by Tree Species",
       x = "Tree Species",
       y = "Count") +
  theme_minimal()
# Sample data frame

##calculating percent survival for each species,

# percent survival calcualtions -------------------------------------------

summary_data_survival <- merged_data_all %>%
  group_by(plot_num, species, innoculation) %>%  
  filter(!is.na(survival)) %>% # Group by both variables
  summarise(
    total = n(),                                   # Total count of rows in the group
    count_Y = sum(survival == "Y"),                # Count of 'Y' in the survival column
    percent_Y = (count_Y / total) * 100,           # Calculate percentage of 'Y'
    .groups = "drop"                               # Ungroup after summarising
  )
summary_data_survival <- summary_data_survival%>%
  group_by(species, innoculation)%>%
  summarise(
    survival = mean(percent_Y),
    sd_survival = sd(percent_Y),
    n = n()
  ) %>%
  mutate(
    se_survival = sd_survival / sqrt(n)
  )


# survival plot -----------------------------------------------------------


ggplot(summary_data_survival, aes(x = species, y = survival, fill = innoculation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +  # Adjust dodge and bar width
  geom_errorbar(
    aes(ymin = survival- se_survival, ymax = survival + se_survival),
    position = position_dodge(width = 0.8),  # Ensure error bars align with bars
    width = 0.25  # Width of the error bars
  ) +
  labs(
    #title = "",
    x = "Species",
    y = "survival",
    fill = "Inoculation Status"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


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

#######

