### Two-stage meta-analysis on pollinator abundance in relation to distance_m from natural habitat

# Load the libraries
library(ggplot2)
library(dplyr)
library(MASS)
library(metafor)
library(here)

# Load the data
setwd(rprojroot::find_rstudio_root_file())
df <- read.csv("raw_data/Fulldataset_abundance.csv")
head(df)
colnames(df)

# Remove rows with NA in 'abundance' or 'distance_m' columns
df <- df[!is.na(df$abundance_all) & !is.na(df$distance_m), ]
nrow(df)

# Get all individual report names
unique_report <- unique(df$report)
length(unique_report) #number of studies
unique_report

# Count unique locations per author and add total sum
author_location_count <- df %>%
  group_by(authors) %>%
  summarise(unique_locations = n_distinct(location)) %>%
  bind_rows(tibble(authors = "Total", unique_locations = sum(.$unique_locations)))
print(author_location_count)

range(df$distance_m, na.rm = TRUE)

# Calculate the range of distance_m values for each report
distance_range_by_report <- df %>%
  group_by(report) %>%
  summarise(min_distance = min(distance_m, na.rm = TRUE),
            max_distance = max(distance_m, na.rm = TRUE),
            n_observations = n()) # Optional: count number of observations per report

# Print the results
print(distance_range_by_report)

# Calculate summary statistics for max_distance
summary_stats <- distance_range_by_report %>%
  summarise(mean_max_distance = mean(max_distance, na.rm = TRUE),
            median_max_distance = median(max_distance, na.rm = TRUE),
            min_max_distance = min(max_distance, na.rm = TRUE),
            max_max_distance = max(max_distance, na.rm = TRUE))

# Print summary statistics
print(summary_stats)

# Merge max_distance into df
df <- df %>%
  left_join(distance_range_by_report %>% dplyr::select(report, max_distance), by = "report")

# View as a sorted table for better readability
distance_range_by_report %>% arrange(min_distance, max_distance)
write.csv(distance_range_by_report, "outputs/abundance/Distance range by report - abundance.csv", row.names = FALSE)


########################## Estimate Effect Sizes ###########################

# Loop through each report and create a separate data frame for each
for (report in unique_report) {
  
  # Filter the df for the current report
  report_data <- df %>% filter(report == !!report)
  
  # Create a variable name dynamically based on the report name
  dataset_name <- paste0("A_", gsub(" ", "_", report))
  
  # Assign the filtered data frame to a new variable in the global environment
  assign(dataset_name, report_data, envir = .GlobalEnv)
}

############################## Fit GLM Models ###############################
# Loop through each report to fit the GLM and store the results
for (report in unique_report) {
  dataset_name <- paste0("A_", report)   # Construct the dataset name
  data <- get(dataset_name)
  
  # Determine model type
  balanced_effort <- length(unique(data$sampling_effort)) > 1
  
  # Print report-level details
  print(paste("Report:", report))
  print(paste("Balanced sampling effort:", balanced_effort))
  
  # Check if sampling effort is constant within the report
  if (balanced_effort) {
    model <- glm.nb(abundance_all ~ log(distance_m + 1) + offset(log(sampling_effort)), data = data, link = log) # Negative Binomial model with offset
  } else {
    model <- glm.nb(abundance_all ~ log(distance_m + 1), data = data, link = log) # Simple Negative Binomial model
  }
  
  # Save the model
  model_name <- paste0("GLM_A_", report)
  assign(model_name, model, envir = .GlobalEnv)
  
  # Print summary for debugging
  print(paste("Model stored as", model_name))
  print(summary(model))
}

# Create an empty data frame to store the results
results <- data.frame(
  Report = character(),
  Authors = character(),
  Slope = numeric(),
  StdError = numeric(),
  PValue = numeric(),
  AgrIntensity = character(),
  Sites = character(),
  Habitat = character(),
  Pollinator = character(),
  Method = character(),
  DistanceMeasure = character(),
  MaxDistance = numeric(),
  RoB = character(),
  stringsAsFactors = FALSE
)

# Loop through each report to extract coefficients and store them
for (report in unique_report) {
  # Construct model name
  model_name <- paste0("GLM_A_", report)
  dataset_name <- paste0("A_", report)
  
  # Get the model and dataset from the global environment
  model <- get(model_name)
  data <- get(dataset_name)
  
  # Extract coefficients and their statistics from the model
  coef_summary <- summary(model)$coefficients
  
  # Extract the slope, standard error, and p-value for the log(distance_m + 1) term
  slope <- coef_summary["log(distance_m + 1)", "Estimate"]
  std_error <- coef_summary["log(distance_m + 1)", "Std. Error"]
  p_value <- coef_summary["log(distance_m + 1)", "Pr(>|z|)"]
  
  # Extract additional information from the dataset
  authors <- unique(data$authors)
  agr_intensity <- unique(data$agr_intensity)
  sites <- unique(data$sites)
  habitat <- unique(data$habitat)
  pollinator <- unique(data$pollinator)
  sampling_method <- unique(data$sampling_method)
  distance_measure <- unique(data$distance_measure)
  max_distance <- unique(data$max_distance)
  #RoB <- unique(data$RoB)
  
  # Append the results to the data frame
  results <- rbind(results, data.frame(
    Authors = authors,
    Slope = slope,
    StdError = std_error,
    PValue = p_value,
    AgrIntensity = agr_intensity,
    Sites = sites,
    Habitat = habitat,
    Pollinator = pollinator,
    Method = sampling_method,
    DistanceMeasure = distance_measure,
    MaxDistance = max_distance
    #RoB = RoB
  ))
}

# Print the overview of the extracted results
print(results)

# Add Variance as a column 
results$Variance <- results$StdError^2

# Save the results to a CSV file --> For Stage 2 of the meta-analysis
write.csv(results, "outputs/abundance/Abundance-distance stage 1 results.csv", row.names = FALSE)

### Plot model fit for each individual report ####

# Define the folder path to save the plots
dir.create(here("outputs", "abundance", "GLM fits"), recursive = TRUE, showWarnings = FALSE)

# Loop through each report to generate plots
for (report in unique_report) {
  # Construct dataset and model names
  dataset_name <- paste0("A_", report)
  model_name <- paste0("GLM_A_", report)
  
  # Get dataset and model
  dataset <- get(dataset_name)
  model <- get(model_name)
  
  # Create a smooth sequence of distances for prediction
  smooth_distance <- data.frame(distance_m = seq(min(dataset$distance_m), max(dataset$distance_m), length.out = 100))
  
  # Check if sampling_effort was used as an offset
  balanced_effort <- "sampling_effort" %in% colnames(dataset) && length(unique(dataset$sampling_effort)) > 1
  
  # If sampling_effort was used as an offset, include it in smooth_distance
  if (balanced_effort) {  
    smooth_distance$sampling_effort <- mean(dataset$sampling_effort, na.rm = TRUE)  # Use a representative value
  }
  
  # Predict values for smooth distance sequence
  preds <- predict(model, newdata = smooth_distance, type = "response", se.fit = TRUE)
  
  # Add predictions and confidence intervals to smooth dataset
  smooth_distance <- smooth_distance %>%
    mutate(
      predicted = preds$fit,
      lower_CI = preds$fit - 1.96 * preds$se.fit,
      upper_CI = preds$fit + 1.96 * preds$se.fit
    )
  
  # Plot raw data with a **smooth** fitted line and confidence intervals
  p <- ggplot(dataset, aes(x = distance_m, y = abundance_all)) +  # Base plot with actual data
    geom_ribbon(data = smooth_distance, aes(x = distance_m, ymin = lower_CI, ymax = upper_CI), 
                inherit.aes = FALSE, fill = "grey70", alpha = 0.4) +  # Confidence interval shading
    geom_line(data = smooth_distance, aes(x = distance_m, y = predicted), 
              inherit.aes = FALSE, color = "blue", linewidth = 1) +  # Smooth fitted line
    geom_point(size = 3, alpha = 0.6, colour = "black") + # Raw data points
    labs(x = "distance in m", y = "Pollinator abundance (all)", title = paste(report, "et al.")) +
    theme_minimal(base_size = 12)  # Clean theme
  
  # Print the plots
  print(p)
  
  # Save the plots
  ggsave(filename = here("outputs", "abundance", "GLM fits", paste0("Model_Fit_", report, ".png")), 
         plot = p, width = 8, height = 6, dpi = 300)}

########################## Meta-analysis ###########################

### Conduct the meta-analysis using metafor
# Load the CSV file into a dataframe
abundance_es <- read.csv("outputs/abundance/Abundance-distance stage 1 results.csv", header=TRUE, stringsAsFactors = FALSE)
colnames(abundance_es)

### Fit a random-effects model using the calculated effect sizes ###
res <- rma(yi = Slope, vi = Variance, data = abundance_es, method = "REML")
res # report tau^2, tau, I^2 and test for heterogeinity (Q test)
predict(res, digits=3)
confint(res) # Look at confidence intervals (useful to report this)

weights(res)

# Create the forest plot
# Set standard font size for metafor
par(font = 1, cex = 1.2)  # Use same scale as ggplot

tiff(filename = here("outputs", "abundance", "Abundance (all) forest plot.tiff"), width = 9, height = 7, units = "in", res = 600)  # Open TIFF device
forest(res,
       slab = abundance_es$Authors,                # Labels for the studies
       xlab = "Slope",                       # Label for the x-axis
       xlim = c(-2.3, 2),                            # Customize x-axis limits
       refline = 0,                                # Add reference line at 0
       header = "a) Total pollinator abundance (all species)",         # Header for the plot
       annotate = TRUE,                            # Add study annotations
       ilab.xpos = -0.015,                         # Adjust position of study effect size labels
       cex = 0.8)                                  # Manage overall font size
dev.off()  # Close device to save the file

############### Plot decay curve ####################

# Get the minimum and maximum distances from the dataset
min_distance <- min(df$distance_m, na.rm = TRUE)
distance <- min_distance:max_distance # Define distances

# Print the min and max distance values
print(paste("Min distance:", min_distance))
print(paste("Max distance:", max_distance))

# Set relative intercept a = 1
a <- 1

# Define slopes from meta-analysis
mean_slope <- res$beta[1]
lower_CI <- res$ci.lb      
upper_CI <- res$ci.ub      

# Create a dataframe for the mean effect and confidence intervals
df_abundance <- data.frame(
  distance = distance,
  rel_abundance_mean = a * exp(mean_slope * log(distance + 1)),  # Mean effect
  rel_abundance_upper = a * exp(upper_CI * log(distance + 1)),  # Upper CI
  rel_abundance_lower = a * exp(lower_CI * log(distance + 1))  # Lower CI
)

# Decay curve
decaycurve <- ggplot(df_abundance, aes(x = distance)) +
  geom_ribbon(aes(ymin = rel_abundance_lower, ymax = rel_abundance_upper), fill = "lightgrey", alpha = 0.5) +  # Shaded CI region
  geom_line(aes(y = rel_abundance_mean), color = "blue", size = 1) +  # Mean effect
  geom_line(aes(y = rel_abundance_upper), linetype = "dashed", color = "black") +  # Upper CI
  geom_line(aes(y = rel_abundance_lower), linetype = "dashed", color = "black") +  # Lower CI
  scale_x_continuous(limits = c(0, 5000), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +  # Force y-axis to start at 0
  labs(
    x = "Distance to nearest natural habitat (m)",
    y = "Relative change in pollinator abundance",
    title = "c) Predicted decay curve") +
  theme_classic(base_size = 14) +  # Standard font size
  theme(
    plot.title = element_text(size = 9, face = "bold"),  # Match forest plot title
    axis.title = element_text(size = 9),  # Match axis labels
    axis.text = element_text(size = 9),  # Standardize tick labels
    plot.margin = margin(10, 20, 10, 10)  # Adjust spacing
  )
decaycurve
ggsave(filename = here("outputs", "abundance", "Abundance decay curve.png"), plot = decaycurve, width = 6, height = 5, units = "in", dpi = 300)

#### Calculate % decline at 1km from natural habitat
# Define the function for predicted abundance
predict_abundance <- function(distance_m, slope) {
  a * exp(slope * log(distance_m + 1))
}

# Calculate predicted abundances at 0m and 1000m
predicted_abundance_0m <- predict_abundance(0, mean_slope)
predicted_abundance_1000m <- predict_abundance(1000, mean_slope)

# Compute percentage decline
percentage_decline <- (1 - (predicted_abundance_1000m / predicted_abundance_0m)) * 100

# Print result
print(paste("Predicted % decline at 1km:", round(percentage_decline, 2), "%"))

########################## Additional analyses ###########################

### Moderator analysis ###
# Summarise each moderators distribution
moderators <- c("Method", "DistanceMeasure" , "MaxDistance", "AgrIntensity", "Pollinator", "Habitat") # List of moderator variables

# Use lapply to get counts for each moderator
moderator_summaries <- lapply(moderators, function(moderator) {
  abundance_es %>%
    count(!!sym(moderator), name = "Report_Count") %>%
    arrange(desc(Report_Count))
})

names(moderator_summaries) <- moderators # Name the list elements
moderator_summaries # Print all summaries

## Meta-analysis with individual moderators
# Moderator analysis for agricultural intensity
res.modintensity <- rma(Slope, Variance, mods = ~ 0 + AgrIntensity, data=abundance_es)
res.modintensity

# Moderator analysis for distance measure
res.modpollinator <- rma(Slope, Variance, mods = ~ 0 + Pollinator, data=abundance_es)
res.modpollinator

# Moderator analysis for habitat type
res.modhabitat <- rma(Slope, Variance, mods = ~ 0 + Habitat, data=abundance_es)
res.modhabitat

########################## Sensitivity analyses ###########################

# Test for small study bias incl publication bias
metafor::funnel(res) #should be symmetric
regtest(res) # Egger's regression test. If this is not statistically significant, no evidence for small study bias
ranktest(res)

### Influential studies ###
inf <- influence(res) # Check if any individual studies were very influential
print(inf) # influential studies have an * next to them
plot(inf) # red dots for influential studies 

# Re-run meta-analysis excluding each study one by one
leave1out <- leave1out(res, digits = 3)
leave1out
write.csv(leave1out, "outputs/abundance/Leave1out abundance.csv", row.names = FALSE)

# QQ-normal plots
qqnorm(res, main = "Random-Effects Model")

### Sensitivity analyses --> methods
res.modmethod <- rma(Slope, Variance, mods = ~ 0 + Method, data=abundance_es)
res.modmethod

res.moddistance_method <- rma(Slope, Variance, mods = ~ 0 + DistanceMeasure, data=abundance_es)
res.moddistance_method

# Add categories for the maximum distance scales (small, medium, large) for sensitivty analysis
abundance_es <- abundance_es %>%
  mutate(DistanceCategory = case_when(
    MaxDistance < 750 ~ "small",
    MaxDistance >= 750 & MaxDistance <= 3000 ~ "medium",
    MaxDistance > 3000 ~ "large"
  ))

# Convert DistanceCategory to a factor
abundance_es$DistanceCategory <- factor(abundance_es$DistanceCategory, levels = c("small", "medium", "large"))

res.modmaxdistance <- rma(Slope, Variance, mods = ~ 0 + DistanceCategory, data=abundance_es)
res.modmaxdistance

### Subgroup analysis excl studies with high RoB ###
# Filter for medium RoB studies
# results_mediumRoB <- subset(abundance_es, RoB == "medium")

# Calculate variance for the meta-analysis
# results_mediumRoB$Variance <- results_mediumRoB$StdError^2

# Fit random-effects model
# res_abundance_mediumRoB <- rma(yi = Slope, vi = Variance, data = results_mediumRoB)
# print(res_abundance_mediumRoB)

# Create a forest plot
# forest(res_abundance_mediumRoB,
#       slab = results_mediumRoB$Authors,                # Labels for the studies
#       xlab = "Slope",                                  # Label for the x-axis
#       xlim = c(-4, 3),                                 # Customize x-axis limits
#       refline = 0,                                     # Add reference line at 0
#       header = "Pollinator abundance subgroup; medium Risk of Bias", # Header for the plot
#       annotate = TRUE,                                 # Add study annotations
#       ilab.xpos = -0.015,                              # Adjust position of study effect size labels
#       cex = 0.8                                        # Manage overall font size
#       )

########################## Wild pollinators ###########################

# Subgroup analysis for wild pollinators only (dataset excl Apis mellifera and Apis cerana where it was managed)
wild_data <- df[!is.na(df$abundance_wild), ]# Filter rows with non-NA values in abundance_wild
nrow(wild_data)

# Get all individual report names
unique_report_wild <- unique(wild_data$report)
unique_report_wild
length(unique_report_wild) #number of studies

# Loop through each report in the wild dataset
for (report in unique_report_wild) {
  # Filter the dataset for the current report
  dataset_name <- paste0("AWild_", report)
  data <- wild_data[wild_data$report == report, ]
  
  # Assign the filtered dataset to a new variable in the global environment
  assign(dataset_name, data, envir = .GlobalEnv)
  
  # Fit the GLM for abundance_wild
  model <- glm.nb(abundance_wild ~ log(distance_m + 1), data = data, link = log)
  
  # Dynamically name the model
  model_name <- paste0("GLM_AWild_", report)
  assign(model_name, model, envir = .GlobalEnv)
  
  # Print model summary
  print(paste("Model stored as", model_name))
  print(summary(model))
}

# Create an empty data frame for results
results_wild <- data.frame(
  Authors = character(),
  Slope = numeric(),
  StdError = numeric(),
  PValue = numeric(),
  AgrIntensity = character(),
  Sites = character(),
  Habitat = character(),
  Pollinator = character(),
  Method = character(),
  DistanceMeasure = character(),
  MaxDistance = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each report in the wild dataset
for (report in unique_report_wild) {
  # Construct model and dataset names
  model_name <- paste0("GLM_AWild_", report)
  dataset_name <- paste0("AWild_", report)
  
  # Get the model and dataset
  model <- get(model_name)
  data <- get(dataset_name)
  
  # Extract coefficients
  coef_summary <- summary(model)$coefficients
  slope <- coef_summary["log(distance_m + 1)", "Estimate"]
  std_error <- coef_summary["log(distance_m + 1)", "Std. Error"]
  p_value <- coef_summary["log(distance_m + 1)", "Pr(>|z|)"]
  
  # Additional information from the dataset
  authors <- unique(data$authors)
  agr_intensity <- unique(data$agr_intensity)
  sites <- unique(data$sites)
  habitat <- unique(data$habitat)
  pollinator <- unique(data$pollinator)
  sampling_method <- unique(data$sampling_method)
  distance_m_measure <- unique(data$distance_m_measure)
  max_distance_m <- unique(data$max_distance_m)
  
  # Append results
  results_wild <- rbind(results_wild, data.frame(
    Authors = authors,
    Slope = slope,
    StdError = std_error,
    PValue = p_value,
    AgrIntensity = agr_intensity,
    Sites = sites,
    Habitat = habitat,
    Pollinator = pollinator,
    Method = sampling_method,
    DistanceMeasure = distance_measure,
    MaxDistance = max_distance
  ))
}

# Save the results
write.csv(results_wild, "outputs/abundance/Wild abundance-distance stage 1 results.csv", row.names = FALSE)

# Calculate variance for the meta-analysis
results_wild$Variance <- results_wild$StdError^2

# Fit random-effects model
res_abundance_wild <- rma(yi = Slope, vi = Variance, data = results_wild)
print(res_abundance_wild)

# Create a forest plot
tiff(file = here("outputs", "abundance", "Abundance (wild) forest plot.tiff"), width = 9, height = 7, units = "in", res = 600)  # Open TIFF device
forest(res_abundance_wild,
       slab = results_wild$Authors,                # Labels for the studies
       xlab = "Slope",                       # Label for the x-axis
       xlim = c(-3, 2),                            # Customize x-axis limits
       refline = 0,                                # Add reference line at 0
       header = "b) Wild pollinator abundance (excluding managed honeybees)",         # Header for the plot
       annotate = TRUE,                            # Add study annotations
       ilab.xpos = -0.015,                         # Adjust position of study effect size labels
       cex = 0.8                                   # Manage overall font size
)    
dev.off()  # Close device to save the file