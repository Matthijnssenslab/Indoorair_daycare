#lineplot of average pcr value for enteroviruses and counts of how many times enteroviruses were detected.

library(dplyr)
library(ggplot2)
library(lubridate)
library(readr) 

# load pcr data
PCR_data <- read_csv(file ="all_air_test_data.csv")

# Convert samplingDate column to date format if it's not already
PCR_data$samplingDate <- as.Date(PCR_data$samplingDate, format = "%Y-%m-%d")

# Filter for human bocavirus and dates within 2022
filtered_pcr <- PCR_data %>%
  filter(pathogen == "human enterovirus (incl. rhinovirus)", year(samplingDate) == 2022, detected== "TRUE")
# Assuming the 'filtered_pcr' data frame is already created

# Group by month and calculate mean, standard deviation, and sample size
monthly_stats <- filtered_pcr %>%
  group_by(month = floor_date(samplingDate, "month")) %>%
  summarise(
    average_Ct = mean(Ct_value, na.rm = TRUE),
    sd_Ct = sd(Ct_value, na.rm = TRUE), # Standard deviation
    n = n() # Sample size
  ) %>%
  mutate(
    SEM = sd_Ct / sqrt(n), # Calculate Standard Error of the Mean (SEM)
    CI_lower = average_Ct - qt(0.975, n - 1) * SEM, # Lower bound of the 95% CI
    CI_upper = average_Ct + qt(0.975, n - 1) * SEM # Upper bound of the 95% CI
  )

# Display the monthly statistics with confidence intervals
print(monthly_stats)


# Plot the average Ct_value per month with confidence intervals
ggplot(monthly_stats, aes(x = month, y = average_Ct)) +
  geom_line(color = "darkblue") +  # Draw line connecting average Ct_values
  geom_point(color = "darkblue") +  # Add points for each month's average Ct_value
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 1, color = "darkblue") +  # Add confidence interval error bars
  labs(title = "Average Ct Values of Detected Human Enterovirus in 2022",
       x = "Month",
       y = "Average Ct Value") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +  # Format x-axis dates
  scale_y_reverse(limits = c(40, 28))  # Invert y-axis and set limits


library(ggplot2)
library(dplyr)
library(lubridate)
library(readr)

######IQRange####
# Group by month and calculate mean and standard deviation
monthly_stats <- filtered_pcr %>%
  group_by(month = floor_date(samplingDate, "month")) %>%
  summarise(
    average_Ct = mean(Ct_value, na.rm = TRUE),
    sd_Ct = sd(Ct_value, na.rm = TRUE), # Standard deviation
    n = n() # Sample size
  ) %>%
  mutate(
    sd_lower = average_Ct - sd_Ct, # Lower bound of the standard deviation
    sd_upper = average_Ct + sd_Ct  # Upper bound of the standard deviation
  )


# Plot the average Ct_value per month with standard deviation error bars
ggplot(monthly_stats, aes(x = month, y = average_Ct)) +
  geom_line(color = "darkblue") +  # Draw line connecting average Ct_values
  geom_point(color = "darkblue") +  # Add points for each month's average Ct_value
  geom_errorbar(aes(ymin = sd_lower, ymax = sd_upper), width = 0.2, color = "gray") +  # Add standard deviation error bars
  labs(title = "Average Ct Values of Detected Human Enterovirus in 2022 with SD",
       x = "Month",
       y = "Average Ct Value") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +  # Format x-axis dates
  scale_y_reverse(limits = c(40, 30))  # Invert y-axis and set limits

print(monthly_stats)



# Group by month and calculate mean, standard deviation, and the mid 50% range
monthly_stats <- filtered_pcr %>%
  group_by(month = floor_date(samplingDate, "month")) %>%
  summarise(
    average_Ct = mean(Ct_value, na.rm = TRUE),
    sd_Ct = sd(Ct_value, na.rm = TRUE), # Standard deviation
    n = n(), # Sample size
    IQR_lower = quantile(Ct_value, 0.25, na.rm = TRUE), # 25th percentile
    IQR_upper = quantile(Ct_value, 0.75, na.rm = TRUE)  # 75th percentile
  )

# Display the monthly statistics with the mid 50% range
print(monthly_stats)

# Plot the average Ct_value per month with the mid 50% range
ggplot(monthly_stats, aes(x = month, y = average_Ct)) +
  geom_line(color = "darkblue") +  # Draw line connecting average Ct_values
  geom_point(color = "darkblue") +  # Add points for each month's average Ct_value
  geom_errorbar(aes(ymin = IQR_lower, ymax = IQR_upper), width = 2, color = "#7E6AAA") +  # Add IQR range bars
  labs(title = "Average Ct Values of Detected Human Enterovirus in 2022 with Mid 50% Range",
       x = "Month",
       y = "Average Ct Value") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +  # Format x-axis dates
  scale_y_reverse(
    limits = c(36, 31), # Set the limits for y-axis
    breaks = seq(31, 36, length.out = 5) # Create 5 breaks between the limits
  )



###counts####

# Filter to include only rows where 'detected' is "TRUE", then group by month and count
monthly_detected_count <- filtered_pcr %>%
  filter(detected == "TRUE") %>% # Filter for samples with 'detected' as "TRUE"
  group_by(month = floor_date(samplingDate, "month")) %>%
  summarise(detected_count = n()) # Count the number of TRUE detections

# Display the monthly detected count
print(monthly_detected_count)
# Group by month and detection status, then count the number of occurrences
monthly_detection_summary <- filtered_pcr %>%
  group_by(month = floor_date(samplingDate, "month"), detection_status = detected) %>%
  summarise(count = n(), .groups = 'drop') # Count the number of occurrences for each detection status

# Display the monthly detection summary
print(monthly_detection_summary)

