#heatmap of Fig 2A

# Load necessary libraries
library(dplyr)
library(tidyr)
library(pheatmap)
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
summary_readcounts <- filtered_data_families  %>%
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


# Filter the data for specifically rhinovirus subspecies
filtered_data_rhinovirus <- filtered_data_families %>%
  dplyr::filter(grepl("rhinovirus", subspecies, ignore.case = TRUE))

# Define all possible months to ensure coverage
all_months <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")

# Aggregate total reads_aligned for each Month and subspecies
aggregated_reads_monthly <- filtered_data_rhinovirus %>%
  group_by(Month, subspecies) %>%
  summarize(total_reads_aligned = sum(reads_aligned, na.rm = TRUE), .groups = 'drop')

# Count distinct sample IDs for each Month across all subspecies
overall_sample_counts <- filtered_data_families %>%
  group_by(Month) %>%
  summarize(overall_sample_count = n_distinct(sample_ID), .groups = 'drop')

# Join aggregated monthly reads with overall sample counts
aggregated_data <- aggregated_reads_monthly %>%
  left_join(overall_sample_counts, by = "Month")

# Calculate reads aligned per sample for each Month and subspecies
aggregated_data <- aggregated_data %>%
  mutate(normalized_reads_per_sample = total_reads_aligned)

# Apply log transformation to normalized_reads_per_sample
aggregated_data$log_normalized_reads_per_sample <- log10(aggregated_data$normalized_reads_per_sample + 1) # Adding 1 to avoid log(0)

# Pivot wider and fill any missing values with 0
heatmap_data_prep <- aggregated_data %>%
  select(subspecies, Month, log_normalized_reads_per_sample) %>%
  pivot_wider(names_from = Month, values_from = log_normalized_reads_per_sample, values_fill = list(log_normalized_reads_per_sample = 0))

# Add missing months with 0 values
for (month in all_months) {
  if (!month %in% names(heatmap_data_prep)) {
    heatmap_data_prep[[month]] <- 0
  }
}

# Ensure the Month columns are in the correct order
heatmap_data_prep <- heatmap_data_prep %>%
  select(subspecies, all_months)

# Convert the data frame to a numeric matrix
heatmap_matrix <- as.matrix(heatmap_data_prep[-1])  # Remove the subspecies column for the matrix
rownames(heatmap_matrix) <- heatmap_data_prep$subspecies

# Plot the heatmap using pheatmap
pheatmap(heatmap_matrix,
         main = "Log Transformed Normalized Reads Aligned per Sample for Rhinovirus Subspecies by Month",
         fontsize = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#70243A"))(100))
# Count distinct sample IDs for each Month where rhinovirus was detected
rhinovirus_detected_counts <- filtered_data_rhinovirus %>%
  group_by(Month) %>%
  summarize(detected_sample_count = n_distinct(sample_ID), .groups = 'drop')

# Print the counts
print(rhinovirus_detected_counts)
# Count distinct sample IDs for each Month across all subspecies (overall sample count)
overall_sample_counts <- filtered_data_families %>%
  group_by(Month) %>%
  summarize(overall_sample_count = n_distinct(sample_ID), .groups = 'drop')

# Join the rhinovirus detected counts with overall sample counts
monthly_sample_counts <- rhinovirus_detected_counts %>%
  left_join(overall_sample_counts, by = "Month")

# Print the counts
print(monthly_sample_counts)
