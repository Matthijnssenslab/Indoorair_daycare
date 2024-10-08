#Figure 5.A

# Load required libraries
library(dplyr)
library(tidyr)  
library(ggplot2)
library(RColorBrewer)
library(readr)

Mastertable <- read.csv2('../codes_to_submit/csv_files/denovo/mastertable_filtered.csv')
# Drop unwanted columns
Mastertable <- Mastertable %>%
  select(-CE17, -CE40, -CE41, -CE53)

# Manually create the vector of sample_IDs in the desired order
sample_order <- c(
  "CE18", "CE31", "CE50", "CE51", "CE52", "CE49", "CE19", "CE11", "CE38",
  "CE37", "CE10", "CE32", "CE34", "CE36", "CE35", "CE30", "CE33", "CE29",
  "CE42", "CE43", "CE44", "CE45", "CE46", "CE47", "CE20", "CE48", "CE28",
  "CE26", "CE15", "CE21", "CE7", "CE14", "CE6", "CE4", "CE22", "CE3", "CE2",
  "CE25", "CE1", "CE39", "CE5", "CE24", "CE13", "CE27", "CE9", "CE12", "CE16",
  "CE8", "CE23"
)
# Normalize counts to get proportions within each sample
superkingdom_proportions <- Mastertable %>%
  select(starts_with("CE"), Final_superkingdom) %>%
  group_by(Final_superkingdom) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(-Final_superkingdom, names_to = "Sample", values_to = "Reads") %>%
  group_by(Sample) %>%
  mutate(Proportion = Reads / sum(Reads)) %>%
  ungroup() %>%
  mutate(Sample = factor(Sample, levels = sample_order)) %>%  # Use factor to set levels according to sample_order
  arrange(Final_superkingdom, Sample)  # Arrange primarily by Superkingdom and then by Sample

# Plot
ggplot(superkingdom_proportions, aes(x = Sample, y = Proportion, fill = Final_superkingdom)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Reads Mapped to Superkingdoms",
       x = "Sample",
       y = "Proportion",
       fill = "Final Superkingdom") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Adjust Final_superkingdom to merge "Other" with "Archaea" and "Root" with "Unknown"
Mastertable3 <- Mastertable %>%
  mutate(Final_superkingdom = case_when(
    Final_superkingdom %in% c("Other", "Archaea") ~ "Other",
    Final_superkingdom %in% c("Root", "Unknown", NA) ~ "Unknown",
    Final_superkingdom %in% c("Unclassified viruses", "Eukaryotic viruses", "Prokaryotic viruses", "Viruses") ~ "Viruses",
    TRUE ~ as.character(Final_superkingdom)
  ))

# Sum up reads for each (adjusted) superkingdom across all samples and calculate proportions
overall_proportions <- Mastertable3 %>%
  select(starts_with("CE"), Final_superkingdom) %>%
  group_by(Final_superkingdom) %>%
  summarise(Total_Reads = sum(across(starts_with("CE")))) %>%
  ungroup() %>%
  mutate(Proportion = Total_Reads / sum(Total_Reads))
# Calculate percentages
overall_proportions$Percentage <- overall_proportions$Proportion * 100

# Plot
ggplot(overall_proportions, aes(x = "", y = Proportion, fill = Final_superkingdom)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = paste0(Total_Reads, " (", round(Percentage, 1), "%)")), 
            position = position_stack(vjust = 0.5)) +
  labs(title = "Overall Proportion of Reads Mapped to Superkingdoms",
       x = "",
       y = "Proportion",
       fill = "Final Superkingdom") +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


# Adjust Final_superkingdom to rename "Viruses" to "Prokaryotic viruses"
Mastertable4 <- Mastertable %>%
  mutate(Final_superkingdom = ifelse(Final_superkingdom == "Viruses", "Prokaryotic viruses", Final_superkingdom))

# Filter for virus-related rows
virus_data <- Mastertable %>%
  filter(Final_superkingdom %in% c("Unclassified viruses", "Eukaryotic viruses", "Prokaryotic viruses"))

# Sum up reads for each virus category across all samples and calculate proportions
virus_proportions <- virus_data %>%
  select(starts_with("CE"), Final_superkingdom) %>%
  group_by(Final_superkingdom) %>%
  summarise(Total_Reads = sum(across(starts_with("CE")))) %>%
  ungroup() %>%
  mutate(Proportion = Total_Reads / sum(Total_Reads), 
         Label = paste0("(", Total_Reads, ", ", scales::percent(Proportion), ")"))

# Plot
ggplot(virus_proportions, aes(x = "", y = Proportion, fill = Final_superkingdom)) +
  geom_bar(stat = "identity", width = 1, show.legend = FALSE) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Overall Proportion of Viral Reads by Category",
       x = "",
       y = "Proportion") +
  scale_fill_brewer(palette = "Dark2") +  # Change the color palette here
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


# Filter for data related to phylums
family_data <- virus_data %>%
  filter(!is.na(Final_family) & Final_superkingdom == "Eukaryotic viruses") # Ensure rows have a phylum defined

# Sum up reads for each family across all samples and calculate proportions
family_proportions <- family_data %>%
  select(starts_with("CE"), Final_family) %>%
  group_by(Final_family) %>%
  summarise(Total_Reads = sum(across(starts_with("CE"))), .groups = 'drop') %>%
  mutate(Proportion = Total_Reads / sum(Total_Reads)) %>%
  # Group families with less than 2% proportion into "Other"
  mutate(Final_family = ifelse((Proportion < 0.02 | Final_family == "Retroviridae") & Final_family != "Papillomaviridae", "Other", Final_family)) %>%
  group_by(Final_family) %>%
  summarise(Total_Reads = sum(Total_Reads), .groups = 'drop') %>%
  mutate(Proportion = Total_Reads / sum(Total_Reads),
         Label = paste0("(", Total_Reads, ", ", scales::percent(Proportion), ")")) %>%
  # Order by Proportion, ensuring "Other" is at the end
  mutate(Final_family = fct_reorder(fct_infreq(Final_family), Proportion, .desc = TRUE, .fun = sum),
         Final_family = fct_relevel(Final_family, "Other", after = Inf))

# Plot
ggplot(family_proportions, aes(x = "", y = Proportion, fill = Final_family)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Overall Proportion of Reads by Family",
       x = "",
       y = "Proportion") +
  scale_fill_brewer(palette = "Set3") +  # Color-blind safe palette
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank())





