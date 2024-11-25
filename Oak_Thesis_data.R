library(readxl)
library(tidyverse)
colonization_data <- read_excel("/Users/zacharyshortt/Desktop/DNA_work_SMM/Colonization_data.xlsx")
master_data <-read_excel("/Volumes/f004qd8/Samo_Data/master_data.xlsx")

combined_coldata <- colonization_data %>%
  group_by(`Tree code`) %>%
  summarize(across(everything(), sum, na.rm = TRUE), .groups = 'drop')


print(combined_coldata)

# Merging the data frames
merged_data_all <- merge(master_data,combined_coldata , by = "Tree code", all = TRUE)

# View the result
print(combined_coldata)


colnames(merged_data_all)

# rename the column since it had a bad name 
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

# Step 2: Calculate radius and volume (cone model)
merged_data_all <- merged_data_all %>%
  mutate(
    radius = width_7_5_24 / 2,  # Radius in cm (width is already in cm)
    volume = (1/3) * pi * (radius^2) * Heightin_6_8_24  # Volume calculation
  )

# Step 4: Plot the volume against species using ggplot2
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

# Step 2: Plot with side-by-side bars and error bars
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



# summary_data_col <- merged_data_all %>%
#   group_by(species, innoculation) %>%
#   summarise(
#     col_avg = mean(percent_col, na.rm = TRUE),  # Mean colonization
#      sd_col= sd(percent_col, na.rm = TRUE),      # Standard deviation
#     n = n()                                    # Number of observations
#   ) %>%
#   mutate(
#     se_col = sd_col / sqrt(n-1)  # Calculate standard error (SE)
#   )

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
#   geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +  # Adjust dodge and bar width
#   geom_errorbar(
#     aes(ymin = col_avg - se_col, ymax = col_avg + se_col),
#     position = position_dodge(width = 0.8),  # Ensure error bars align with bars
#     width = 0.25  # Width of the error bars
#   ) +
#   labs(
#     title = "Tree colonization by Species and Inoculation Status with Error Bars",
#     x = "Species",
#     y = "Mean Colonization (%)",
#     fill = "Inoculation Status"
#   ) +
#   theme_minimal() +
#   theme(legend.position = "top")

ggplot(summary_data_col_plot, aes(x = species, y = col_avg, fill = plot_num)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_errorbar(
    aes(ymin = col_avg - se_col, ymax = col_avg + se_col),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  labs(
    title = "Tree colonization by Species and Inoculation Status with Error Bars",
    x = "Species",
    y = "Mean Colonization (%)",
    fill = "Plot Number"  # Adjust legend title
  ) +
  scale_fill_manual(values = c("1b" = "grey", "2b" = "grey", "2a" = "grey", "3"= "blue", "4" = "red", "5" = "red", "6" = "blue", "7"= "red", "8" = "blue")) +  # Customize colors
  theme_minimal() +
  theme(legend.position = "top")

# Filter out plot numbers 1b, 2a, and 2b
summary_data_col_plot_filtered <- summary_data_col_plot %>%
  filter(!plot_num %in% c("1b", "2a", "2b"))

# Create the plot
ggplot(summary_data_col_plot_filtered, aes(x = species, y = col_avg, fill = plot_num)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_errorbar(
    aes(ymin = col_avg - se_col, ymax = col_avg + se_col),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  labs(
    title = "Tree colonization by Species and Inoculation Status with Error Bars",
    x = "Species",
    y = "Mean Colonization (%)",
    fill = "Plot Number"  # Adjust legend title
  ) +
  scale_fill_manual(values = c("3"= "blue", "4" = "red", "5" = "red", "6" = "blue", "7"= "red", "8" = "blue")) +  # Customize colors
  theme_minimal() +
  theme(legend.position = "top")





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

# 
# ##Calcualting volume of tree 
# 
# # Step 1: Calculate Volume as a new column
# merged_data_all <- merged_data_all %>%
#   mutate(
#     radius = width_7_5_25 / 2,  # Calculate radius
#     volume = (1/3) * pi * (radius^2) * Heightin_6_8_24  # Calculate volume
#   )
# 
# # Step 2: Display the first few rows to check the new column
# head(merged_data_all)
# 
# # Step 3: Plot the volume against species using ggplot2
# ggplot(merged_data_all, aes(x = species, y = volume, fill = species)) +
#   geom_bar(stat = "identity") +
#   labs(
#     title = "Tree Volume by Species",
#     x = "Species",
#     y = "Volume (Cone Model)"
#   ) +
#   theme_minimal() +
#   theme(legend.position = "none")
# 
