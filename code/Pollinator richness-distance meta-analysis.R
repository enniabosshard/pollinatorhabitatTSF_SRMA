### Two-stage meta-analysis on pollinator richness in relation to isolation from natural habitat

# Load the libraries
library(ggplot2)
library(dplyr)
library(MASS)
library(lme4)
library(metafor)
library(ggeffects)
library(here)

# Load data
setwd(rprojroot::find_rstudio_root_file())
df <- read.csv("raw_data/Fulldataset_richness.csv", header = TRUE, stringsAsFactors = FALSE)

# Remove rows with NA in 'richness' or 'distance' columns
df <- df[!is.na(df$richness_all) & !is.na(df$distance_m), ]
nrow(df)
colnames(df)

# Get all individual report names
unique_report <- unique(df$report)
length(unique_report) #number of reports
unique_report

# Check distribution of study design types across datasets
df %>%
  distinct(report, study_design) %>%  # keep only unique study-design pairs
  count(study_design)  # count how many studies per design

# Count unique locations (i.e. study sites) per author and add total sum
author_location_count <- df %>%
  group_by(authors) %>%
  summarise(unique_locations = n_distinct(location)) %>%
  bind_rows(tibble(authors = "Total", unique_locations = sum(.$unique_locations)))
print(author_location_count)

# Summarise species richness per study
 richness_overview <- df %>%
  group_by(study) %>%
  summarise(
    min_species = min(richness_all, na.rm = TRUE),
    max_species = max(richness_all, na.rm = TRUE),
    mean_species = mean(richness_all, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(richness_overview, "outputs/richness/Study-level species richness overview.csv", row.names = FALSE)

range(df$distance_m, na.rm = TRUE)

# Calculate the range of distance values for each report
distance_range_by_report <- df %>%
  group_by(report) %>%
  summarise(min_distance = min(distance_m, na.rm = TRUE),
            max_distance = max(distance_m, na.rm = TRUE),
            n_observations = n()) 

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
write.csv(distance_range_by_report, "outputs/richness/Study-level distance ranges.csv", row.names = FALSE)

########################## Estimate Effect Sizes ###########################

## Estimate effect sizes for each individual report
# Loop through each report and create a separate data frame for each
for (report in unique_report) {
  
  # Filter the df for the current report
  report_data <- df %>% filter(report == !!report)
  
  # Create a variable name dynamically based on the report name
  dataset_name <- paste0("R_", gsub(" ", "_", report))
  
  
  # Assign the filtered data frame to a new variable in the global environment
  assign(dataset_name, report_data, envir = .GlobalEnv)
}

## Fit Models ###
# Loop through each report to fit the models and store the results
for (report in unique_report) {
  # Construct the dataset name
  dataset_name <- paste0("R_", report)
  data <- get(dataset_name)
  
  # Determine model type
  unbalanced_effort <- length(unique(data$repeat_measures)) > 1
  design <- unique(data$study_design)
  
  # Print report-level details
  print(paste("Study:", report))
  print(paste("Study design:", design))
  print(paste("Unbalanced sampling effort:", unbalanced_effort))
  
  # Choose model type based on study design and sampling effort
  if (design == "single distance per site") {
    # Use GLM
    if (unbalanced_effort) {
      model <- glm.nb(richness_all ~ log(distance_m + 1) + offset(log(repeat_measures)), data = data, link = log)
    } else {
      model <- glm.nb(richness_all ~ log(distance_m + 1), data = data, link = log)
    }
    
  } else {
    # Use GLMM with random effect for location
    if (unbalanced_effort) {
      model <- glmer.nb(richness_all ~ log(distance_m + 1) + offset(log(repeat_measures)) + (1 | location), data = data)
    } else {
      model <- glmer.nb(richness_all ~ log(distance_m + 1) + (1 | location), data = data)
    }
  }
  
  # Save the model
  model_name <- paste0("model_R_", report)
  assign(model_name, model, envir = .GlobalEnv)
  
  # Print summary for debugging
  print(paste("Model stored as", model_name))
  print(summary(model))
}

# Create an empty data frame to store the results
results <- data.frame(
  Authors = character(),
  Slope = numeric(),
  StdError = numeric(),
  PValue = numeric(),
  AgrIntensity = character(),
  Sites = character(),
  Design = character(),
  Habitat = character(),
  Pollinator = character(),
  Method = character(),
  DistanceMeasure = character(),
  MaxDistance = numeric(),
  TaxonomicResolution = character(),
  stringsAsFactors = FALSE
)

# Loop through each report to extract coefficients and store them
for (report in unique_report) {
  # Construct model name
  model_name <- paste0("model_R_", report)
  dataset_name <- paste0("R_", report)
  
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
  design <- unique(data$study_design)
  habitat <- unique(data$habitat)
  pollinator <- unique(data$pollinator)
  sampling_method <- unique(data$sampling_method)
  distance_measure <- unique(data$distance_measure)
  max_distance <- unique(data$max_distance)
  taxonomic_resolution <- unique(data$taxonomic_resolution)
  
  # Append the results to the data frame
  results <- rbind(results, data.frame(
    Authors = authors,
    Slope = slope,
    StdError = std_error,
    PValue = p_value,
    AgrIntensity = agr_intensity,
    Sites = sites,
    Design = design,
    Habitat = habitat,
    Pollinator = pollinator,
    Method = sampling_method,
    DistanceMeasure = distance_measure,
    MaxDistance = max_distance,
    TaxonomicResolution = taxonomic_resolution
  ))
}

# Print the overview of the extracted results
print(results)

# Save the results to a CSV file --> For Stage 2 of the meta-analysis
write.csv(results, "outputs/richness/Richness-distance stage 1 results.csv", row.names = FALSE)

### Plot model fit for each individual report ####
# Create folder to save the plots
dir.create(here("outputs", "richness", "model fits"), recursive = TRUE, showWarnings = FALSE)

# Loop through each report to generate plots
for (report in unique_report) {
  # Construct dataset and model names
  dataset_name <- paste0("R_", report)
  model_name <- paste0("model_R_", report)
  
  # Get dataset and model
  dataset <- get(dataset_name)
  model <- get(model_name)
  
  # Create smooth sequence of distances
  smooth_distance <- data.frame(distance_m = seq(min(dataset$distance_m), max(dataset$distance_m), length.out = 100))
  
  ### For the GLMs
  if (inherits(model, "glm")) {

    # Check if repeat_measures used as offset
    balanced_effort <- "repeat_measures" %in% colnames(dataset) && length(unique(dataset$repeat_measures)) > 1
    if (balanced_effort) {
      smooth_distance$repeat_measures <- mean(dataset$repeat_measures, na.rm = TRUE)
    }
    
    # Predict with SEs
    preds <- predict(model, newdata = smooth_distance, type = "response", se.fit = TRUE)
    smooth_distance <- smooth_distance %>%
      mutate(
        predicted = preds$fit,
        lower_CI = preds$fit - 1.96 * preds$se.fit,
        upper_CI = preds$fit + 1.96 * preds$se.fit
      )
    
    # Build plot
    p <- ggplot(dataset, aes(x = distance_m, y = richness_all)) +
      geom_ribbon(data = smooth_distance, aes(x = distance_m, ymin = lower_CI, ymax = upper_CI),
                  inherit.aes = FALSE, fill = "grey70", alpha = 0.4) +
      geom_line(data = smooth_distance, aes(x = distance_m, y = predicted),
                inherit.aes = FALSE, color = "blue", linewidth = 1) +
      geom_point(size = 3, alpha = 0.6, colour = "black") +
      labs(x = "Distance to natural habitat (m)",
           y = "Pollinator richness (all)",
           title = paste(report, "et al. (GLM)")) +
      theme_minimal(base_size = 12)
    
    ggsave(filename = here("outputs", "richness", "model fits",
                           paste0(report, "_GLM_fit.png")),
           plot = p, width = 8, height = 6, dpi = 300)
    
    next
  }
  
  ### For the GLMMs
  if (inherits(model, "glmerMod")) {
    
    # Predictions from ggeffects package
    preds <- ggeffects::ggemmeans(model, terms = list(distance_m = smooth_distance$distance_m), 
                                  condition = c(repeat_measures = mean(dataset$repeat_measures)))  
    
    # Build plot (same style)
    p <- ggplot() +
      geom_ribbon(data = preds, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey70", alpha = 0.4) +
      geom_line(data = preds, aes(x = x, y = predicted), colour = "blue", linewidth = 1) +
      geom_point(data = dataset, aes(x = distance_m, y = richness_all), colour = "black", alpha = 0.6, size = 3) +
      labs(x = "Distance to natural habitat (m)", y = "Pollinator richness (all)", title = paste(report, "et al. (GLMM)")) +
      theme_minimal(base_size = 12)
    
    ggsave(filename = here("outputs", "richness", "model fits",
                           paste0(report, "_GLMM_fit.png")),
           plot = p, width = 8, height = 6, dpi = 300)
  }
}

########################## Meta-Analysis ###########################

### Conduct the meta-analysis using metafor
richness_es <- results
richness_es$Variance <- richness_es$StdError^2

### Fit a random-effects model using the calculated effect sizes ###
res <- rma(yi = Slope, vi = Variance, data = richness_es)
res # report tau^2, tau, I^2 and test for heterogeinity (Q test)
predict(res, digits=3)
confint(res) # Look at confidence intervals (useful to report this)

# Create the forest plot
tiff(filename = here("outputs", "richness", "Richness (all) forest plot.tiff"), width = 9, height = 7, units = "in", res = 600)  # Open TIFF device
forest(res,
       slab = richness_es$Authors,                    # Labels for the reports
       xlab = "Slope",               # Label for the x-axis
       xlim = c(-2, 2),                      # Customize x-axis limits if needed
       refline = 0,                                # Add reference line at 0
       header = "a) Total pollinator richness (all species)",                           # Header for the plot
       annotate = TRUE,                            # Add report annotations
       ilab.xpos = -0.015,                         # Adjust position of report effect size labels
       cex = 0.8)                                  # Manage overall font size
dev.off()  # Close device to save the file

##### Plot decay surve
# Define distances (log scale for consistency with Ricketts et al.)

# Get the minimum and maximum distances from the dataset
min_distance <- min(df$distance_m, na.rm = TRUE)
# Create a sequence of distances for plotting (e.g., 100 points between min and max)
distance_seq <- seq(min_distance, max_distance, length.out = 100)

# Print the min and max distance values
print(paste("Min distance:", min_distance))
print(paste("Max distance:", max_distance))

# Set relative intercept r = 1
r <- 1

# Define slopes from meta-analysis
mean_slope <- res$beta[1]
lower_CI <- res$ci.lb      
upper_CI <- res$ci.ub

# Create a dataframe for the mean effect and confidence intervals
df_richness <- data.frame(
  distance = distance_seq,
  richness_mean = r * exp(mean_slope * log(distance_seq + 1)),  # Mean effect
  richness_upper = r * exp(upper_CI * log(distance_seq + 1)),  # Upper CI
  richness_lower = r * exp(lower_CI * log(distance_seq + 1))  # Lower CI
)

# Plot using ggplot2
decaycurve <- ggplot(df_richness, aes(x = distance_seq)) +
  geom_ribbon(aes(ymin = richness_lower, ymax = richness_upper), fill = "lightgrey", alpha = 0.5) +  # Shaded CI region
  geom_line(aes(y = richness_mean), color = "blue", size = 1) +  # Mean effect
  geom_line(aes(y = richness_upper), linetype = "dashed", color = "black") +  # Upper CI
  geom_line(aes(y = richness_lower), linetype = "dashed", color = "black") +  # Lower CI
  scale_x_continuous(breaks = seq(0, 5000, by = 1000), limits = c(0, 5000), expand = c(0, 0)) +  # Force x-axis to start at 0
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +  # Force y-axis to start at 0
  labs(
    x = "Distance to nearest natural habitat (m)",
    y = "Relative change in pollinator richness",
    title = "b) Predicted decay curve") +
  theme_classic(base_size = 14) +  # Standard font size
  theme(
    plot.title = element_text(size = 9, face = "bold"),  # Match forest plot title
    axis.title = element_text(size = 9),  # Match axis labels
    axis.text = element_text(size = 9),  # Standardize tick labels
    plot.margin = margin(10, 20, 10, 10)  # Adjust spacing
  )
decaycurve
ggsave(filename = here("outputs", "richness", "Richness decay curve.png"), plot = decaycurve, width = 6, height = 5, units = "in", dpi = 300)

#### Calculate % decline at 1km from natural habitat
# Define the function for predicted richness
predict_richness <- function(distance, slope) {
  r * exp(slope * log(distance + 1))
}

# Calculate predicted richnesss at 0m and 1000m
predicted_richness_0m <- predict_richness(0, mean_slope)
predicted_richness_1000m <- predict_richness(1000, mean_slope)

# Compute percentage decline
percentage_decline <- (1 - (predicted_richness_1000m / predicted_richness_0m)) * 100

# Print result
print(paste("Predicted % decline at 1km:", round(percentage_decline, 2), "%"))

########################## Additional analyses ###########################

### Moderator analysis ###
# Summarise each moderators distribution
moderators <- c("Method", "DistanceMeasure" , "AgrIntensity", "Pollinator", "Habitat", "TaxonomicResolution") # List of moderator variables

# Use lapply to get counts for each moderator
moderator_summaries <- lapply(moderators, function(moderator) {
  richness_es %>%
    count(!!sym(moderator), name = "Report_Count") %>%
    arrange(desc(Report_Count))
})

names(moderator_summaries) <- moderators # Name the list elements
moderator_summaries # Print all summaries

## Meta-analysis with individual moderators
# Moderator analysis for habitat type
res.modhabitat <- rma(Slope, Variance, mods = ~ 0 + Habitat, data=richness_es)
res.modhabitat

# Moderator analysis for agricultural intensity
res.modintensity <- rma(Slope, Variance, mods = ~ 0 + AgrIntensity, data=richness_es)
res.modintensity

########################## Sensitivity analyses ###########################

# Test for small study bias incl publication bias
metafor::funnel(res) #should be symmetric
regtest(res) # Egger's regression test. If this is not statistically significant, no evidence for small study bias
ranktest(res)

### Influential reports ###
inf <- influence(res) # Check if any individual reports were very influential
print(inf) # influential reports have an * next to them
plot(inf) # red dots for influential reports 

# Re-run meta-analysis excluding each report one by one
leave1out <- leave1out(res, digits = 3)
leave1out
write.csv(leave1out, ("outputs/richness/Leave1out_richness.csv"), row.names = FALSE)

# QQ-normal plots
qqnorm(res, main = "Random-Effects Model")

## Sensitivity analyses
# Sensitivity model 1: Pollinator sampling method
res.modmethod <- rma(Slope, Variance, mods = ~ 0 + Method, data=richness_es)
res.modmethod

# Sensitivty model 2: Distance measure
res.moddistance <- rma(Slope, Variance, mods = ~ 0 + DistanceMeasure, data=richness_es)
res.moddistance

# Sensitivity model 3: Distance category (max distance)
# Add categories for the maximum distance scales (small, medium, large) for sensitivty analysis
richness_es <- richness_es %>%
  mutate(DistanceCategory = case_when(
    MaxDistance < 750 ~ "small",
    MaxDistance >= 750 & MaxDistance <= 3000 ~ "medium",
    MaxDistance > 3000 ~ "large"
  ))

# Convert DistanceCategory to a factor
richness_es$DistanceCategory <- factor(richness_es$DistanceCategory, levels = c("small", "medium", "large"))

res.modmaxdistance <- rma(Slope, Variance, mods = ~ 0 + DistanceCategory, data=richness_es)
res.modmaxdistance

# Sensitivity model 4: Taxonomic resolution
res.modtaxonomic <- rma(Slope, Variance, mods = ~ 0 + TaxonomicResolution, data=richness_es)
res.modtaxonomic
