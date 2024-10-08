#Figure 4.A


library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(readr)

filtered_viruses <- read.csv("filtered_data_animal.csv")

# Your provided order of sample IDs
sample_order <- c("CE18", "CE31", "CE50", "CE51", "CE52", "CE49", "CE11", "CE38", "CE37", "CE10", "CE32", 
                  "CE34", "CE36", "CE35", "CE30", "CE33", "CE29", "CE43", "CE44", "CE46", "CE47", "CE28", 
                  "CE26", "CE15", "CE21", "CE7", "CE14", "CE6", "CE4", "CE22", "CE3", "CE2", "CE25", "CE39", 
                  "CE5", "CE24", "CE13", "CE27", "CE9", "CE12", "CE16", "CE23", "CE17", "CE53", "CE40", "CE41")
# Aggregate certain species under the 'Parvoviridae' family
filtered_viruses <- filtered_viruses %>%
  mutate(
    aggregate_family = case_when(
      species %in% c("Avian dependoparvovirus 1", "Galliform chaphamaparvovirus 2", "Galliform aveparvovirus 1") ~ "Parvoviridae",
      TRUE ~ as.character(species)
    )
  )

# Ensure all combinations of aggregate_family and sample_ID are included
complete_data <- expand.grid(aggregate_family = unique(filtered_viruses$aggregate_family), sample_ID = sample_order)

# Merge with the filtered viruses to include read counts, filling missing values with 0
merged_data <- merge(complete_data, filtered_viruses, by = c("aggregate_family", "sample_ID"), all.x = TRUE)
merged_data$reads_aligned[is.na(merged_data$reads_aligned)] <- 0  # Replace NA with 0

# Reshape the data to create the abundance matrix
abundance_matrix <- dcast(merged_data, aggregate_family ~ sample_ID, value.var = "reads_aligned", fun.aggregate = sum)

# Apply log10 transformation to the abundance values (adding 1 to avoid log(0))
abundance_matrix[,-1] <- log10(abundance_matrix[,-1] + 1)

# Prepare the matrix for the heatmap
aggregate_names <- abundance_matrix$aggregate_family
heatmap_matrix <- as.matrix(abundance_matrix[,-1])
row.names(heatmap_matrix) <- aggregate_names

# Reorder columns according to the provided sample order
heatmap_matrix <- heatmap_matrix[, sample_order]

# Define a color palette for the heatmap
color_palette <- colorRampPalette(c( "white", "#84ACAC"))(n = 100)

# Generate the heatmap
pheatmap(heatmap_matrix,
         color = color_palette,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "none")
