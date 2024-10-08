#Human associated viruses, de novo assembly.
#steps

# Required libraries
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# data <- read.csv("~/Mastertable_filtered.csv")
# Filter Mastertable (filtered for contaminants) for "Eukaryotic viruses"
Eukaryotic_viruses <- Mastertable %>%
  filter(Final_superkingdom == "Eukaryotic viruses")

#associated genera and families were selected to focus on
selected_genera <- c("Unclassified Polyomaviridae", "Alphapolyomavirus quintihominis", "Betapapillomavirus", "Betapolyomavirus","Gammapapillomavirus","Bocaparvovirus","Cytomegalovirus", "Deltapolyomavirus", "Dependoparvovirus", "Picobirnavirus", "Unclassified Papillomaviridae","Alphapolyomavirus"  )

selected_families <- c("Anelloviridae")


# Filtering steps
Eukaryotic_viruses_filtered <- Eukaryotic_viruses %>%
  filter(Final_genus %in% selected_genera | Final_family %in% selected_families) 

# Update Final_genus based on merging criteria
Eukaryotic_viruses_filtered <- Eukaryotic_viruses_filtered %>%
  mutate(Final_genus = case_when(
    Final_genus %in% c("Gammapapillomavirus", "Unclassified Papillomaviridae") ~ "Gammapapillomavirus",
    Final_genus %in% c("Deltapolyomavirus", "Unclassified Polyomaviridae") ~ "Deltapolyomavirus",
    Final_genus %in% c("Unclassified Anelloviridae", "Betatorquevirus") ~ "Anelloviridae",
    TRUE ~ Final_genus
  ))

# Ensure only existing columns in sample_order are used
valid_samples <- intersect(sample_order, colnames(Eukaryotic_viruses_filtered))

# Aggregate the data by 'Final_genus' and 'Group', summing up the valid samples
aggregated_data <- Eukaryotic_viruses_filtered %>%
  group_by(Final_genus) %>%
  summarise(across(all_of(valid_samples), sum, na.rm = TRUE)) %>%
  ungroup()

# Order the aggregated data by 'Group' and then by 'Final_genus'
aggregated_data <- aggregated_data %>%
  arrange(Final_genus)


# Create the abundance matrix and the row annotation for 'Group'
abundance_matrix <- as.matrix(aggregated_data %>% select(all_of(valid_samples)))
rownames(abundance_matrix) <- aggregated_data$Final_genus

# Convert the abundance values to numeric
abundance_matrix <- apply(abundance_matrix, 2, as.numeric)

# Log-transform the abundance data
abundance_data_log <- log10(abundance_matrix + 1)

# Create row annotations for the heatmap
row_annotation <- aggregated_data %>%
  select(Final_genus) %>%
  distinct() %>%
  arrange( Final_genus)

# Ensure row order in the heatmap matches the row annotation
#abundance_data_log_or <- abundance_data_log[match(row_annotation$Final_genus, rownames(abundance_data_log)), ]
abundance_data_log <- abundance_data_log[complete.cases(abundance_data_log), ]
# Set the row names of abundance_data_log to 'Final_genus' values
rownames(abundance_data_log) <- row_annotation$Final_genus
# Define a custom color palette for log-scaled values
custom_color_palette <- colorRampPalette(c("gray99", "pink", "pink4"))(100)

# Set the row names of abundance_data_log to 'Final_genus' values
rownames(abundance_data_log) <- row_annotation$Final_genus

# Plot the heatmap with ordered rows, 'Group' annotation, and custom colors
pheatmap(abundance_data_log, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         annotation_row = row_annotation,
         color = custom_color_palette)

# Define color breaks for the legend
legend_breaks = seq(-2, 2, length.out = 101) 

# Generate labels for each break
legend_labels = round(seq(min(legend_breaks), max(legend_breaks), length.out = length(legend_breaks)), 1)

# Match the order of the samples with the column names of the abundance matrix
col_annotation <- metadata %>%
  select(Names, Month) %>%
  filter(Names %in% colnames(abundance_data_log)) %>%
  arrange(match(Names, colnames(abundance_data_log)))

# Ensure the row names of the annotation are the same as the column names of the abundance data
rownames(col_annotation) <- col_annotation$Names
# Convert the 'Month' column to a factor
col_annotation$Month <- factor(col_annotation$Month, levels = unique(col_annotation$Month))

# Generate a color palette that matches the number of unique months
# Ensure that the number of colors matches the number of unique months
month_colors <- setNames(brewer.pal(n = length(unique(col_annotation$Month)), name = "Set3"), unique(col_annotation$Month))

# Now create the heatmap with the corrected annotation_colors
pheatmap(
  abundance_data_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  color = custom_color_palette,
  main = "Heatmap Title",
  fontsize_row = 8,
  fontsize_col = 8,
  legend = TRUE,
  legend_breaks = legend_breaks,
  legend_labels = legend_labels,
)

#Combining figures and further edits were done using Adobe Illustrator.