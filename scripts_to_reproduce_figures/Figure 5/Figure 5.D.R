#Figure 5.D

library(phyloseq)
library(readr)
library(tidyr)
library(dplyr)


filtered_airfilters <- readRDS("input/filtered_airfilters.rds")

# Load metadata
metadata <- read_csv("csv_files/denovo/metadata2.csv")


metadata <- metadata %>%
  # Select necessary columns
  select(Names, Month, Sample_or_Control) %>%
  
  # Filter out controls
  filter(Sample_or_Control != "Control") %>%
  
  # Define the two groups based on months
  mutate(Group = ifelse(Month %in% c("April", "May", "June", "July", "August", "September", "October"), 
                        "April-October", "Other Months"))

library(dplyr)
library(ggplot2)
# Melt the phyloseq object into a data frame
family_counts <- psmelt(filtered_airfilters) %>%
  filter(Final_superkingdom == "Eukaryotic viruses" & 
           !Final_family %in% c("Retroviridae", "Phycodnaviridae", "Unclassified Cressdnaviricota")) %>%
  group_by(Sample, Final_family) %>%
  summarize(Family_Count = sum(Abundance), .groups = 'drop')

# Joining family_counts with metadata
full_data <- family_counts %>%
  rename(Sample = Sample) %>%  # Rename if necessary
  inner_join(metadata, by = c("Sample" = "Names"))

# Filter out specified families
full_data <- full_data %>%
  filter(!Final_family %in% c("Parvoviridae", "Polyomaviridae", "Anelloviridae", "Other", "Papillomaviridae","Picobirnaviridae","Herpesviridae","Mimiviridae","Circoviridae"))


# Filter to include only rows where Family_Count is greater than 0
filtered_data <- full_data %>%
  filter(Family_Count > 0)



# Aggregate data to count the number of unique samples where each family was detected by group
detections_count <- filtered_data %>%
  group_by(Final_family, Group) %>%
  summarise(detection_count = n_distinct(Sample), .groups = 'drop') %>%
  arrange(Final_family, Group)

# Summarize the total detection counts for each group
grouped_detections_count <- detections_count %>%
  group_by(Group) %>%
  summarise(total_detections = sum(detection_count), .groups = 'drop',
            total_undetections = 130-sum(detection_count))

# Print the detection counts for each group
print(grouped_detections_count)
# Perform the chi-squared test
chi_square_test <- chisq.test(grouped_detections_count[, c("total_detections", "total_undetections")])

# Print the results
print(chi_square_test)


# Print the results with a more readable p-value
p_value <- format.pval(chi_square_test$p.value, digits = 4)
print(paste("Chi-Square Test p-value:", p_value))

# Define significance level
significance_level <- 0.05
significance_label <- ifelse(chi_square_test$p.value < significance_level, "***", "ns")
# Find the maximum y-value for positioning the annotation
max_y <- max(grouped_detections_count$total_detections)

# Define a custom palette with 18 distinct colors
color_palette <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#B15928"
)
# Plotting the detection counts for the two groups with distinct colors for each family
plot4 <- ggplot(detections_count, aes(x = Group, y = detection_count, fill = Final_family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(title = "Total Detections per Family by Group",
       x = "Group",
       y = "Detection Count",
       fill = "Family") +
  scale_fill_manual(values = color_palette, name = "Family", 
                    breaks = unique(detections_count$Final_family),  # Set breaks to ensure all families appear
                    labels = unique(detections_count$Final_family)) +  # Set labels to match breaks
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 10, face = "bold"),
        panel.grid.major = element_line(color = "black", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 2))+  # Adjust the number of columns in legend
  coord_cartesian(ylim = c(0, 90))  # Set y-axis limits# Adjust the number of columns in legend

# Adjust the number of columns in legend

# Display the plot
print(plot4)



# Plotting the detection counts as stacked bar plot
plot2 <- ggplot(detections_count, aes(x = reorder(Final_family, detection_count), y = detection_count, fill = Group)) +
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
print(plot2)

