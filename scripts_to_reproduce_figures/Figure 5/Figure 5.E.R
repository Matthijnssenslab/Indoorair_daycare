# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(reshape2)

# Read the filtered master table and mapped reads
Mastertable_filtered <- read.csv2('../csv_files/denovo/mastertable_filtered.csv')
Mapped_reads <- read.csv2('../csv_files/denovo/mapped_reads.csv')  # Assuming this file exists
metadata <- read_csv("../csv_files/denovo/metadata2.csv")

# Filter the data for Malasseziaceae family
filtered_mycota <- Mastertable_filtered %>% 
  filter(Final_family == "Malasseziaceae")

# Melt the dataframe to long format (sample-wise reads)
long_df <- melt(filtered_mycota, id.vars = c("CONTIGS", "Final_family"), 
                variable.name = "Sample", value.name = "Reads")

# Convert reads to numeric and remove NAs
long_df$Reads <- as.numeric(long_df$Reads)
long_df <- na.omit(long_df)

# Define sample order and ensure all samples are present
ordered_cols <- c("CE18", "CE31", "CE50", "CE51", "CE52", "CE49", "CE11", "CE38",
                  "CE37", "CE10", "CE32", "CE34", "CE36", "CE35", "CE30", "CE33", 
                  "CE29", "CE43", "CE44", "CE46", "CE47", "CE28", "CE26", "CE15", 
                  "CE21", "CE7", "CE14", "CE6", "CE4", "CE22", "CE3", "CE2", 
                  "CE25", "CE39", "CE5", "CE24", "CE13", "CE27", "CE9", "CE12", 
                  "CE16", "CE23", "CE17", "CE53", "CE40", "CE41")

# Join the sample data with the long format dataframe to ensure all samples are included
sample_df <- data.frame(Sample = ordered_cols)
long_df <- right_join(long_df, sample_df, by = "Sample")

# Ensure correct sample ordering
long_df$Sample <- factor(long_df$Sample, levels = ordered_cols)

# Aggregate total reads by sample
aggregated_df <- long_df %>%
  group_by(Sample) %>%
  summarise(Total_Reads = sum(Reads, na.rm = TRUE), .groups = 'drop')

# Join aggregated reads with mapped reads data to calculate RPM
aggregated_df <- left_join(aggregated_df, Mapped_reads, by = "Sample") %>%
  mutate(RPM = (Total_Reads / Mapped_reads) * 1e6)

# Ensure correct sample ordering for plotting
aggregated_df$Sample <- factor(aggregated_df$Sample, levels = ordered_cols)

# Plot
ggplot(aggregated_df, aes(x = Sample, y = RPM)) +
  geom_bar(stat = "identity", fill = "#E1DB8E") +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = "Sample", y = "RPM", title = "Reads per Million Mapped Reads for Malasseziaceae") +
  scale_y_reverse()
