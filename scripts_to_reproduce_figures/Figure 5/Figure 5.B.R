#Figure 5.B

#load libraries
library(phyloseq)
library(ggplot2)
library(dplyr)

filtered_airfilters <- readRDS("/input/filtered_airfilters.rds")

# Melt the phyloseq object into a data frame
family_counts <- psmelt(filtered_airfilters) %>%
  filter(Final_superkingdom == "Eukaryotic viruses" &
           !Final_family %in% c("Retroviridae", "Phycodnaviridae", "Unclassified Cressdnaviricota")) %>%
  group_by(Sample, Final_family) %>%
  summarize(Family_Count = sum(Abundance), .groups = 'drop')



# Get the overall top 10 families by abundance
top_families <- family_counts %>%
  group_by(Final_family) %>%
  summarize(Total_Abundance = sum(Family_Count), .groups = 'drop') %>%
  top_n(10, Total_Abundance) %>%
  pull(Final_family)

# Recategorize families not in the top 10 as "Other"
family_counts <- family_counts %>%
  mutate(Final_family = if_else(Final_family %in% top_families, as.character(Final_family), "Other")) %>%
  group_by(Sample, Final_family) %>%
  summarize(Family_Count = sum(Family_Count), .groups = 'drop') %>%
  mutate(Relative_Abundance = Family_Count / sum(Family_Count))

# Join with the sample_data to include 'Month' information
sample_data_df <- data.frame(sample_data(filtered_airfilters))
family_counts <- left_join(family_counts, sample_data_df, by = c("Sample" = "Names"))
# Assuming 'family_counts' already has 'Relative_Abundance' calculated per sample

# Ensure the Sample column is ordered according to 'ordered_cols'
family_counts$Sample <- factor(family_counts$Sample, levels = ordered_cols)

# Plotting with the specified order
ggplot(family_counts, aes(x = Sample, y = Relative_Abundance, fill = Final_family)) +
  geom_bar(stat = "identity", position = "fill") +  # 'position = "fill"' to show relative abundance
  theme_minimal() +
  labs(title = "Relative Abundance of Virus Families per Sample", x = "Sample", y = "Relative Abundance") +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  # Adjust text angle and alignment for readability









