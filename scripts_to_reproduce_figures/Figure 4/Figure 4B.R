#This script is intended to create parts of Figure 4b. However, actual visual  was produced in Adobe illustrator by combining parts of the first level and second level below

# Required libraries
library(dplyr)
library(readr)  
library(pheatmap)
library(RColorBrewer)

#####load the mastertable here with a data frame name "Mastertable", )######

 Mastertable <- read.csv('input/mastertable_filtered.csv')
  metadata <- read.csv('csv_files/denovo/metadata2.csv')
  data <-read.csv('csv_files/denovo/denovo_contigs_1.csv')
# Step 1: Filter Mastertable for "Eukaryotic viruses"
Eukaryotic_viruses <- Mastertable %>%
  filter(Final_superkingdom == "Eukaryotic viruses")

selected_genera <- c("Unclassified Partitiviridae", "Lolavirus", "Potexvirus", 
                     "Betachrysovirus", "Unclassified Dicistroviridae", 
                     "Alphaendornavirus", "Unclassified Genomoviridae", "Glossinavirus", 
                     "Chloriridovirus", "Unclassified Nanoviridae", "Lambdapapillomavirus", 
                     "Betapartitivirus", "Deltapartitivirus", "Unclassified Partitiviridae",
                     "Waikavirus", "Unclassified Parvoviridae", "Unclassified Tombusviridae",
                     "Machlomovirus", "Victorivirus", "Unclassified Totiviridae", "Tobamovirus", 
                     "Unclassified Picornavirales", "Unclassified Densovirinae", "Unclassified Parvoviridae", "Gammapolyomavirus", "Ambidensovirus","Alphanucleorhabdovirus","Alternavirus")
# Define the desired order of Sample_IDs
sample_order <- c("CE18", "CE31", "CE50", "CE51", "CE52", "CE49", "CE11", "CE38", "CE37", "CE10", "CE32", 
                  "CE34", "CE36", "CE35", "CE30", "CE33", "CE29", "CE43", "CE44", "CE46", "CE47", "CE28", 
                  "CE26", "CE21", "CE7", "CE14", "CE6", "CE4", "CE22", "CE3", "CE2", "CE25", "CE39", 
                  "CE24", "CE13", "CE27", "CE9", "CE12", "CE16", "CE23", "CE17", "CE53", "CE40", "CE41")

# Filtering steps
Eukaryotic_viruses_filtered <- Eukaryotic_viruses %>%
  filter(Final_genus %in% selected_genera) 

Eukaryotic_viruses_filtered <- left_join(Eukaryotic_viruses_filtered, data, by = c("Final_genus" = "Genus"))

# Ensure only existing columns in sample_order are used
valid_samples <- intersect(sample_order, colnames(Eukaryotic_viruses_filtered))

# Aggregate the data by 'Final_genus' and 'Group', summing up the valid samples
aggregated_data <- Eukaryotic_viruses_filtered %>%
  group_by(Final_genus, Group) %>%
  summarise(across(all_of(valid_samples), sum, na.rm = TRUE)) %>%
  ungroup()

# Order the aggregated data by 'Group' and then by 'Final_genus'
aggregated_data <- aggregated_data %>%
  arrange(Group, Final_genus)


# Create the abundance matrix and the row annotation for 'Group'
abundance_matrix <- as.matrix(aggregated_data %>% select(all_of(valid_samples)))
rownames(abundance_matrix) <- aggregated_data$Final_genus

abundance_matrix <- abundance_matrix[complete.cases(abundance_matrix), ]
abundance_matrix <- apply(abundance_matrix, 2, as.numeric)

# Log-transform the abundance data
abundance_data_log <- log10(abundance_matrix + 1)

# Create row annotations for the heatmap
row_annotation <- aggregated_data %>%
  select(Final_genus, Group) %>%
  distinct() %>%
  arrange(Group, Final_genus)

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
    )

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

# Dynamically set breaks based on the data
breaks <- seq(min(abundance_data_log, na.rm = TRUE), 
              max(abundance_data_log, na.rm = TRUE), 
              length.out = 100)

# Remove custom legend breaks and labels, and allow pheatmap to handle it
pheatmap(
  abundance_data_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  color = custom_color_palette, 
  breaks = breaks,            # Use dynamically generated breaks
  main = "Heatmap of Eukaryotic Virus Genera Abundance",
  fontsize_row = 8,
  fontsize_col = 8)


# this level only creates heatmap of viruses infecting other entities (except animals: plant, fungi or other )
# plant,
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Load the data
data <- read.csv("/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/codes_to_submit/csv_files/denovo/denovo_contigs_1.csv")

# Step 1: Filter Mastertable for "Eukaryotic viruses"
Eukaryotic_viruses <- Mastertable %>%
  filter(Final_superkingdom == "Eukaryotic viruses")

# Define the selected families
selected_families <- c("Alphaflexiviridae", "Chrysoviridae", "Endornaviridae", "Genomoviridae", "Partitiviridae", "Secoviridae", "Tombusviridae", "Totiviridae", "Virgaviridae")

# Filtering steps
Eukaryotic_viruses_filtered <- Eukaryotic_viruses %>%
  filter(Final_family %in% selected_families)

# Ensure only existing columns in sample_order are used
valid_samples <- intersect(sample_order, colnames(Eukaryotic_viruses_filtered))

# Aggregate the data by 'Final_family', summing up the valid samples
aggregated_data <- Eukaryotic_viruses_filtered %>%
  group_by(Final_family) %>%
  summarise(across(all_of(valid_samples), sum, na.rm = TRUE)) %>%
  ungroup()

# Order the aggregated data by 'Final_family'
aggregated_data <- aggregated_data %>%
  arrange(Final_family)

# Create the abundance matrix
abundance_matrix <- as.matrix(aggregated_data %>% select(all_of(valid_samples)))
rownames(abundance_matrix) <- aggregated_data$Final_family

# Convert the abundance values to numeric
abundance_matrix <- apply(abundance_matrix, 2, as.numeric)

# Log-transform the abundance data
abundance_data_log <- log10(abundance_matrix + 1)

# Create row annotations for the heatmap
row_annotation <- aggregated_data %>%
  select(Final_family) %>%
  distinct() %>%
  arrange(Final_family)

# Ensure row order in the heatmap matches the row annotation
abundance_data_log <- abundance_data_log[complete.cases(abundance_data_log), ]
rownames(abundance_data_log) <- row_annotation$Final_family

# Define a custom color palette for log-scaled values
custom_color_palette <- colorRampPalette(c("gray99", "pink", "pink4"))(100)

# Set the row names of abundance_data_log to 'Final_family' values
rownames(abundance_data_log) <- row_annotation$Final_family

# Define color breaks for the legend
legend_breaks <- seq(-2, 2, length.out = 101) 
legend_labels <- round(seq(min(legend_breaks), max(legend_breaks), length.out = length(legend_breaks)), 1)

# Generate the column annotation from metadata
col_annotation <- metadata %>%
  select(Names, Month) %>%
  filter(Names %in% colnames(abundance_data_log)) %>%
  arrange(match(Names, colnames(abundance_data_log)))

# Ensure the row names of the annotation are the same as the column names of the abundance data
rownames(col_annotation) <- col_annotation$Names
col_annotation$Month <- factor(col_annotation$Month, levels = unique(col_annotation$Month))

# Generate a color palette that matches the number of unique months
month_colors <- setNames(brewer.pal(n = length(unique(col_annotation$Month)), name = "Set3"), unique(col_annotation$Month))

# Create the heatmap with the corrected annotation colors
pheatmap(
  abundance_data_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  color = custom_color_palette,
  main = "Heatmap of Selected Virus Families",
  fontsize_row = 8,
  fontsize_col = 8,
  legend = TRUE,
  legend_breaks = legend_breaks,
  legend_labels = legend_labels
)

