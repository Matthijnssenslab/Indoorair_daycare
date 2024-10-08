
#librariess
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra) # for arranging multiple plots
library(pheatmap)
library(tidyverse)
library(ComplexHeatmap)
library(circlize) # Dependency of ComplexHeatmap
library(grid)
library(viridis)
library(scales) # Load the scales package



#load tsv file
data_tsv <- read_tsv("filtered_data_families_human.csv", 
                     col_types = cols())

# Load csv metadata file and ensure dates are read as Date objects
data_csv <- read_csv("combined_metadata.csv")

#metadata and 
filtered_data_families  <- left_join(data_tsv, data_csv, by = c("sample_ID" = "Names"))

#summarize for the total read counts
summary_readcounts <- filtered_data_families%>%
  group_by(sample_ID) %>%
  summarize(
    total_filtered_reads_in_sample = first(total_filtered_reads_in_sample)
  )

mean(summary_readcounts$total_filtered_reads_in_sample, na.rm = TRUE)

median(summary_readcounts$total_filtered_reads_in_sample, na.rm = TRUE)
# Manually create the vector of sample_IDs in the desired order
sample_order <- c(
  "CE18", "CE31", "CE50", "CE51", "CE52", "CE49", "CE19", "CE11", "CE38",
  "CE37", "CE10", "CE32", "CE34", "CE36", "CE35", "CE30", "CE33", "CE29",
  "CE42", "CE43", "CE44", "CE45", "CE46", "CE47", "CE20", "CE48", "CE28",
  "CE26", "CE15", "CE21", "CE7", "CE14", "CE6", "CE4", "CE22", "CE3", "CE2",
  "CE25", "CE1", "CE39", "CE5", "CE24", "CE13", "CE27", "CE9", "CE12", "CE16",
  "CE8", "CE23", "CE17", "CE53", "CE40", "CE41"
)

# Make sure the sample_ID in your data frame is a factor
filtered_data_families$sample_ID <- as.factor(filtered_data_families$sample_ID)

# Set the levels of the sample_ID factor in your data frame according to the provided order
filtered_data_families$sample_ID <- factor(filtered_data_families$sample_ID, levels = sample_order)


# Define all possible months to ensure coverage
all_months <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")

# Filter the data for species Rotavirus A
filtered_data_rotavirus_A <- filtered_data_families %>%
  filter(species == "Rotavirus A")

# Summarize total reads_aligned for each Month where Rotavirus A was detected
rotavirus_reads_aligned <- filtered_data_rotavirus_A %>%
  group_by(Month) %>%
  summarize(total_reads_aligned = sum(reads_aligned, na.rm = TRUE), .groups = 'drop')

# Count distinct sample IDs for each Month across all data (overall sample count)
overall_sample_counts <- filtered_data_families %>%
  group_by(Month) %>%
  summarize(overall_sample_count = n_distinct(sample_ID), .groups = 'drop')

# Join the Rotavirus A reads_aligned with overall sample counts
monthly_reads_aligned <- rotavirus_reads_aligned %>%
  left_join(overall_sample_counts, by = "Month")

# Calculate normalized reads_aligned per sample
monthly_reads_aligned <- monthly_reads_aligned %>%
  mutate(normalized_reads_aligned = total_reads_aligned)

# Ensure all months are included and in the correct order
monthly_reads_aligned$Month <- factor(monthly_reads_aligned$Month, levels = all_months)
monthly_reads_aligned_complete <- expand_grid(Month = factor(all_months, levels = all_months)) %>%
  left_join(monthly_reads_aligned, by = "Month") %>%
  replace_na(list(normalized_reads_aligned = 0))

# Bar plot of normalized reads_aligned for Rotavirus A
ggplot(monthly_reads_aligned_complete, aes(x = Month, y = log10(normalized_reads_aligned), fill = Month)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Normalized Reads Aligned of Rotavirus A per Month",
       x = "Month",
       y = "Normalized Reads Aligned") +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

