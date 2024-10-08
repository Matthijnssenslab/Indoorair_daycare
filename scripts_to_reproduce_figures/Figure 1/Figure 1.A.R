#Figure 1.A
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

#load filtered_data_families_human from the reference_guided folder. This file includes all viruses surviving initial filtering of our bioinformatics analysis. However, below you can find further filtering for human viruses.

# Assuming the initial dataframe is named 'filtered_data_families'
# Our list of viruses and infection types (this list was produced manually by looking at all viruses identified and possibly infecting humans, see supplementary file 1)
virus_list <- data.frame(
  Infection_Type = c("Enteric", "Respiratory", "Other","Other","Other", "Enteric", "Enteric", 
                     "Respiratory", "Skin-associated", "Skin-associated", "Skin-associated", 
                     "Skin-associated", "Skin-associated", "Respiratory", "Respiratory", 
                     "Other", "Respiratory", "Respiratory", "Respiratory", "Respiratory", 
                     "Skin-associated","Skin-associated", "Skin-associated", "Skin-associated", "Skin-associated", 
                     "Skin-associated", "Enteric","Respiratory", "Enteric", "Respiratory"),
  Family = c("Adenoviridae", "Adenoviridae", "Anelloviridae", "Anelloviridae","Anelloviridae", "Astroviridae", "Caliciviridae", 
             "Coronaviridae", "Herpesviridae", "Herpesviridae", "Papillomaviridae", 
             "Papillomaviridae", "Papillomaviridae", "Parvoviridae", "Parvoviridae", "Parvoviridae", "Paramyxoviridae", "Picornaviridae", 
             "Picornaviridae", "Pneumoviridae", "Polyomaviridae","Polyomaviridae", "Polyomaviridae", "Polyomaviridae", "Polyomaviridae", "Poxviridae", 
             "Sedoreoviridae", "Coronaviridae", "Astroviridae", "Coronaviridae"),
  Species_Genus = c("Human mastadenovirus F", "Human mastadenovirus C", "Alphatorquevirus","Betatorquevirus", "unclassified Anelloviridae genus", "Human astrovirus 4", 
                    "Sapporo virus", "Human coronavirus HKU1", "Human betaherpesvirus 5", 
                    "Human alphaherpesvirus 3", "Gammapapillomavirus", "Betapapillomavirus", 
                    "Alphapapillomavirus", "Primate bocaparvovirus 1", "Primate bocaparvovirus 2", 
                    "Adeno-associated dependoparvovirus A", "Human respirovirus 3", "Rhinovirus A", 
                    "Rhinovirus C", "Human orthopneumovirus", "Deltapolyomavirus undecihominis","Deltapolyomavirus decihominis", 
                    "Betapolyomavirus quartihominis", "Betapolyomavirus quartihominis", 
                    "Alphapolyomavirus quintihominis", "Molluscum contagiosum virus", 
                    "Rotavirus A", "Betacoronavirus 1", "Mamastrovirus 1", "Betacoronavirus")
)


# Filter the dataframe for the relevant viruses
filtered_viruses <- filtered_data_families %>%
  filter(species %in% virus_list$Species_Genus | genus %in% virus_list$Species_Genus)

# Group similar genera together using regular expressions. Papillomaviruses and unclassified anelloviruses are grouped together on genus or family level.
filtered_viruses <- filtered_viruses %>%
  mutate(
    genus_aggregate = case_when(
      grepl("^Gammapapillomavirus", species) ~ "Gammapapillomavirus",
      grepl("^Betapapillomavirus", species)  ~ "Betapapillomavirus",
      grepl("^Alphapapillomavirus", species) ~ "Alphapapillomavirus",
      grepl("^Anelloviridae", family) ~ "Anelloviridae",
      TRUE ~ as.character(species)
    )
  )

# Identify all genera under Anelloviridae 
anelloviridae_genera <- unique(filtered_data$genus[filtered_data$family == "Anelloviridae"])

# Append these genera to the virus_list data frame
anelloviridae_df <- data.frame(
  Infection_Type = rep(NA, length(anelloviridae_genera)),
  Family = rep("Anelloviridae", length(anelloviridae_genera)),
  Species_Genus = anelloviridae_genera
)

# Combine virus_list and anelloviridae_df
virus_list <- rbind(virus_list, anelloviridae_df)

joined_df <- filtered_viruses %>%
  left_join(virus_list, by = c("genus_aggregate" = "Species_Genus"), 
            suffix = c("_filtered", "_virus_list"), 
            copy = "shared", 
            keep =TRUE, 
            na_matches = "na",
            relationship = "many-to-many")


# Reshape the data to create a matrix suitable for the heatmap
heatmap_data <- dcast(joined_df, Species_Genus ~ sample_ID, value.var = "reads_aligned", fun.aggregate = sum)

# Apply log10 transformation, assuming all values are positive
heatmap_data[,-1] <- log10(heatmap_data[,-1] + 1)

# Define the desired order of Sample_IDs
ordered_cols <- c("CE18", "CE31", "CE50", "CE51", "CE52", "CE49", "CE11", "CE38", "CE37", "CE10", "CE32", 
                  "CE34", "CE36", "CE35", "CE30", "CE33", "CE29", "CE43", "CE44", "CE46", "CE47", "CE28", 
                  "CE26", "CE21", "CE7", "CE14", "CE6", "CE4", "CE22", "CE3", "CE2", "CE25", "CE39", 
                  "CE24", "CE13", "CE27", "CE9", "CE12", "CE16", "CE23", "CE17", "CE53", "CE40", "CE41")

# Remove the Species_Genus column and convert to matrix
species <- heatmap_data$Species_Genus

# Adjusting ordered_cols to only include columns present in heatmap_data
adjusted_ordered_cols <- ordered_cols[ordered_cols %in% colnames(heatmap_data)]

# Now create the heatmap matrix using the adjusted list of column names
heatmap_matrix <- as.matrix(heatmap_data[, adjusted_ordered_cols])

# Set the row names for the heatmap matrix
row.names(heatmap_matrix) <- heatmap_data$Species_Genus

# Extracting Infection_Type for annotation
annotation_df <- unique(joined_df[, c("Species_Genus", "Infection_Type")])
annotation_df <- annotation_df[!duplicated(annotation_df$Species_Genus), ]
# Ensure Species_Genus is a character and remove any leading/trailing spaces
annotation_df$Species_Genus <- trimws(as.character(annotation_df$Species_Genus))

# Convert to lowercase for matching to avoid case sensitivity issues
annotation_df$Species_Genus <- tolower(annotation_df$Species_Genus)

# Ensure row names of heatmap_matrix are characters and trimmed
row.names(heatmap_matrix) <- trimws(as.character(row.names(heatmap_matrix)))

# Also convert to lowercase
row.names(heatmap_matrix) <- tolower(row.names(heatmap_matrix))

# Now try the match function again
row_annotation <- annotation_df[match(row.names(heatmap_matrix), annotation_df$Species_Genus), "Infection_Type"]

row_annotation <- annotation_df[match(row.names(heatmap_matrix), annotation_df$Species_Genus), "Infection_Type"]
row_annotation_data <- data.frame(Infection_Type = row_annotation)
# row.names(row_annotation_data) <- row.names(heatmap_matrix)

# Convert Infection_Type to a factor and order rows based on this factor
row_annotation_data$Infection_Type <- factor(row_annotation_data$Infection_Type, 
                                             levels = unique(row_annotation_data$Infection_Type))
ordered_indices <- order(row_annotation_data$Infection_Type)
# Ensure you have the tools package available for toTitleCase
row.names(heatmap_matrix) <- tools::toTitleCase(row.names(heatmap_matrix))

ordered_heatmap_matrix <- heatmap_matrix[ordered_indices, ]


# Find the maximum value in the matrix
max_value <- max(heatmap_matrix)

# Define a color palette from blue to red
color_palette <- colorRampPalette(c("#FDFEFF", "#0658A8"))(n = 20)

# Now integrate this col_annotation_matrix in your heatmap
pheatmap(ordered_heatmap_matrix, 
         annotation_row = row_annotation_data[ordered_indices, , drop = FALSE],
         show_rownames = TRUE, 
         show_colnames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = color_palette
)


#Combining figures and further edits were done using Adobe Illustrator.