#Figure 3_A-M only. Phylogenetic trees were c
#Coverage plot parameters can be adjusted in original function.
#choose whichever reference you would like to plot as coverage graph. 
library(rtracklayer)
library(ggplot2)
library(GenomicRanges)
library(dplyr)

#change the seqname depending on the virus you would like to make a coverage plot. Ideally you can optimize it in a way that you do not need to change it, I have not done it yet :)
plotCoverage <- function(bed_file_path, seqname_filter = "PP135058.1_Human_respiratory_syncytial_virus_B", width_filter = 100, bin_size = 100, cap_value = 1000) {
  bed_data <- rtracklayer::import(bed_file_path)
  bed_df <- as.data.frame(bed_data)
  
  # Filter the data frame for only the sequence of interest and width greater than width_filter
  filtered_bed <- bed_df %>%
    filter(seqnames == seqname_filter, width > width_filter)
  
  # Convert the filtered BED data to a GRanges object
  gr <- GRanges(
    seqnames = Rle(filtered_bed$seqnames),
    ranges = IRanges(start = filtered_bed$start, end = filtered_bed$end),
    strand = Rle(filtered_bed$strand)
  )
  
  # Calculate coverage
  cov <- coverage(gr)
  
  # Convert the coverage Rle object to a data frame
  cov_df <- as.data.frame(cov)
  # Reset the row names to get the position information
  cov_df$position <- as.numeric(rownames(cov_df))
  rownames(cov_df) <- NULL
  
  # Bin the positions into bin_size nucleotide intervals
  cov_df$bin <- with(cov_df, cut(position, breaks=seq(from=min(position), to=max(position), by=bin_size), include.lowest=TRUE, labels = FALSE))
  
  # Aggregate the depth by bin
  agg_cov <- aggregate(value ~ bin, data=cov_df, FUN=mean)
  
  # Create midpoints for each bin to use as the x-axis in the plot
  agg_cov$midpoint <- (agg_cov$bin - 0.5) * bin_size
  
  # Cap the highest value to cap_value
  agg_cov$value <- pmin(agg_cov$value, cap_value)
  
  # Get the lowest and highest midpoint values
  lowest_midpoint <- min(agg_cov$midpoint)
  highest_midpoint <- max(agg_cov$midpoint)
  
  # Create the plot (limits have to be changed depending on the virus you would like to plot, otherwise you can use lowest and highest midpoint values. However, this will not be useful if you have a partial genome)
  p <- ggplot(agg_cov, aes(x = midpoint, y = value)) +
    geom_area(stat = "identity", fill = "#106612", alpha = 0.90) +  # Use stat = "identity" for bar plots
    labs(x = paste("Genomic position (binned in", bin_size, "bp intervals)"), y = "Depth of coverage") +
    scale_x_continuous(breaks = seq(from=0, to=20000, by = 1000)) +  # Adjust the x-axis breaks as needed
    coord_cartesian(xlim = c(0, 20000)) +  # Set x-axis limits
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 65, hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid.major = element_line(color = "gray", size = 0.2),
          panel.grid.minor = element_line(color = "gray", size = 0.1))
  
  
  return(p)
}

#astrovirus
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE28.coverage.bed"
coverage_plot_astrovirus <- plotCoverage(plot_path)
print(coverage_plot_astrovirus)

#WU-polyomavirus GU296386.1                       
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE21.coverage.bed"
coverage_plot_polyomavirus <- plotCoverage(plot_path)
print(coverage_plot_polyomavirus)

#papillomavirus-HPV151
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE22.coverage.bed"
coverage_plot_papillomavirus <- plotCoverage(plot_path)
print(coverage_plot_papillomavirus)

#merkel cell
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE13.coverage.bed"
coverage_plot_merkel <- plotCoverage(plot_path)
print(coverage_plot_merkel)

#feline papilloma
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE21.coverage.bed"
coverage_plot_feline <- plotCoverage(plot_path)
print(coverage_plot_feline)

#rhinovirus C43 - OK649386.1_Rhinovirus_C43 
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE24.coverage.bed"
coverage_plot_rhino_c <- plotCoverage(plot_path)
print(coverage_plot_rhino_c)

#rhinovirus A - KY629935.1_Rhinovirus_A59  
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE34.coverage.bed"
coverage_plot_rhino_a <- plotCoverage(plot_path)
print(coverage_plot_rhino_a)

#canary polyoma
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE27.coverage.bed"
coverage_plot_canary <- plotCoverage(plot_path)
print(coverage_plot_canary)

#rsv-b - PP135058.1_Human_respiratory_syncytial_virus_B
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE49.coverage.bed"
coverage_plot_rsv <- plotCoverage(plot_path)
print(coverage_plot_rsv)

plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE51.coverage.bed"
coverage_plot_rsv <- plotCoverage(plot_path)
print(coverage_plot_rsv)


#rota - OR771998.1_Rotavirus_A_VP3_gene,_complete_cds 
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE47.coverage.bed"
coverage_plot_rota <- plotCoverage(plot_path)
print(coverage_plot_rota)

#aav-2 MK139298.1_AAV
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE52.coverage.bed"
coverage_plot_aav <- plotCoverage(plot_path)
print(coverage_plot_aav)

#bocavirus
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE44.coverage.bed"
coverage_plot_boca <- plotCoverage(plot_path)
print(coverage_plot_boca)

#rotaanimal - MF768262.1_Rotavirus_G_NSP4

plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE47.coverage.bed"
coverage_plot_rotag <- plotCoverage(plot_path)
print(coverage_plot_rotag)


#creche densovirus
plot_path <- "/Users/mustafakaratas/Library/CloudStorage/OneDrive-KULeuven/BE/LabVM/airsamples_creche/results_42samples_CE/paper_high_coverage/CE46.coverage.bed"
coverage_plot_denso <- plotCoverage(plot_path)
print(coverage_plot_denso)

