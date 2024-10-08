#Figure 5.C - Human associated viruses testing and plotting included.
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

#load file
filtered_data_families <- read_csv('filtered_data_families_human.csv')
# list of viruses and infection types, choose viruses that we plotted, most likely infecting humans
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


# Combined data filtering and grouping
filtered_data <- filtered_data_families %>%
  filter(species %in% virus_list$Species_Genus | genus %in% virus_list$Species_Genus) %>%
  filter(Sample_or_Control != "Control") %>%
  mutate(Group = ifelse(Month %in% c("April", "May", "June", "July", "August", "September", "October"), 
                        "April-October", "Other Months"))

# Count detections based on `reads_aligned`, considering it detected if `reads_aligned` > 0
detections <- filtered_data %>%
  group_by(sample_ID, family, Group) %>%
  summarise(detected = sum(reads_aligned > 0), .groups = 'drop') %>%
  filter(detected > 0) %>%
  group_by(family, Group) %>%
  summarise(detection_count = n_distinct(sample_ID), .groups = 'drop') %>%
  arrange(family, Group)

# Summarize the total detection counts for each group
grouped_detections_count <- detections %>%
  group_by(Group) %>%
  summarise(total_detections = sum(detection_count), .groups = 'drop',
            total_undetections = 141-sum(detection_count))

# Perform the chi-squared test
chi_square_test <- chisq.test(grouped_detections_count[, c("total_detections", "total_undetections")])

# Print the results
significance_label <- print(p_value)

# Print the results with a more readable p-value
p_value <- format.pval(chi_square_test$p.value, digits = 4)
print(paste("Chi-Square Test p-value:", p_value))

# Find the maximum y-value for positioning the annotation
max_y <- max(grouped_detections_count$total_detections)

color_palette <- c(
  "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
  "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6",
  "#ffff99", "#b15928", "#7fcdbb", "#fdb462", "#386cb0",
  "#f0027f", "#bf5b17", "#666666", "#1b9e77", "#d95f02"
)
# Plotting the detection counts for the two groups with distinct colors for each family
plot3 <- ggplot(detections, aes(x = Group, y = detection_count, fill = family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(title = "Total Detections per Family by Group",
       x = "Group",
       y = "Detection Count",
       fill = "Family") +
  scale_fill_manual(values = color_palette) +  # Use the custom color palette
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10, face = "bold"),
        panel.grid.major = element_line(color = "black", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom") +
  annotate("text", x = 1.5, y = max_y + 2, label = significance_label, size = 6, color = "red")+
  guides(fill = guide_legend(ncol = 2))  +  # Adjust the number of columns in legend
  coord_cartesian(ylim = c(0, 90))  # Set y-axis limits# Adjust the number of columns in legend


# Display the plot
print(plot3)

# Plotting the detection counts as stacked bar plot
plot <- ggplot(detections, aes(x = reorder(family, detection_count), y = detection_count, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = detection_count), 
            position = position_stack(vjust = 0.5), 
            size = 3, 
            color = "white", 
            fontface = "bold") +  # Add bold font to labels
  labs(title = "Total Detections per Family by Group",
       x = "Family",
       y = "Detection Count",
       fill = "Group") +
  scale_fill_manual(values = c("April-October" = "#227A7A", "Other Months" = "#933D4D"), name = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 10, face = "bold"),
        panel.grid.major = element_line(color = "black", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom")

# Display the plot
print(plot)


##Combine figures####
# install.packages("patchwork")
library(ggplot2)
library(patchwork)

# Assuming plot and plot2 are already defined as per previous explanations
# Modify plot1 (plot)
# Combine plots side by side with equal width bars
combined_plot <- plot + plot3 +
  plot_layout(ncol = 2)  # Arrange plots in 2 columns

# Display the combined plot
print(combined_plot)



