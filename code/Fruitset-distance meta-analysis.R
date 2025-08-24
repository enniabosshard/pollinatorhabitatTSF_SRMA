### Two-stage meta-analysis on effects of distance from natural habitat on fruit set

# Load the libraries
library(ggplot2)
library(dplyr)
library(MASS)
library(lme4)
library(ggeffects)
library(metafor)
library(here)
library(glmmTMB)

# Load data
setwd(rprojroot::find_rstudio_root_file())
df <- read.csv("raw_data/Fulldataset_fruitset.csv", header = TRUE, stringsAsFactors = FALSE)

# Remove rows with NA in 'fruitset' or 'distance' columns
df <- df[!is.na(df$fruitset) & !is.na(df$distance_m), ]
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

# Count unique locations per author and add total sum
author_location_count <- df %>%
  group_by(authors) %>%
  summarise(unique_locations = n_distinct(location)) %>%
  bind_rows(tibble(authors = "Total", unique_locations = sum(.$unique_locations)))
print(author_location_count)

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
write.csv(distance_range_by_report,  "outputs/fruitset/Study-level distance ranges.csv", row.names = FALSE)

# Check the distribution of fruit set values
hist(df$fruitset)

# Estimating the study-level effect sizes using beta regression cannot deal with 0s and 1s,
# therefore check how many 0s and 1s in the dataset
sum(df$fruitset == 0, na.rm = TRUE)   # number of true 0s = 11
sum(df$fruitset == 1, na.rm = TRUE)   # number of true 1s = 2

# Truncate the exact 0s and 1s
df$fruitset[df$fruitset == 0] <- 0.00001   # tiny >0
df$fruitset[df$fruitset == 1] <- 0.99999   # just <1

########################## Estimate Effect Sizes ###########################

## Estimate effect sizes for each individual report
# Loop through each report and create a separate data frame for each
for (report in unique_report) {
  
  # Filter the df for the current report
  report_data <- df %>% filter(report == !!report)
  
  # Create a variable name dynamically based on the report name
  dataset_name <- paste0("FS_", gsub(" ", "_", report))
  hist(report_data$fruitset)
  
  # Assign the filtered data frame to a new variable in the global environment
  assign(dataset_name, report_data, envir = .GlobalEnv)
}

## Fit Models ###
# Loop through each report to fit the models and store the results
for (report in unique_report) {
  dataset_name <- paste0("FS_", report) # Construct the dataset name
  data <- get(dataset_name)
  
  # Determine model type
  unbalanced_effort <- length(unique(data$repeat_measures)) > 1
  design <- unique(data$study_design)
  
  # Print report-level details
  print(paste("report:", report))
  print(paste("Study design:", design))
  print(paste("Unbalanced sampling effort:", unbalanced_effort))
  
  # Choose model type based on study design and sampling effort
  if (design == "single distance per site") {
    # Beta regression without random effect
    if (unbalanced_effort) {
      model <- glmmTMB(fruitset ~ log(distance_m + 1), weights = repeat_measures, 
                       family = beta_family(link = "logit"), data = data)
    } else {
      model <- glmmTMB(fruitset ~ log(distance_m + 1), family = beta_family(link = "logit"), data = data)
    }
    
  } else {
    # Beta regression with location as random effect
    if (unbalanced_effort) {
      model <- glmmTMB(fruitset ~ log(distance_m + 1) + (1 | location), weights = repeat_measures, 
                       family = beta_family(link = "logit"), data = data)
    } else {
      model <- glmmTMB(fruitset ~ log(distance_m + 1) + (1 | location), family = beta_family(link = "logit"), data = data)
    }
  }
  
  # Save the models
  model_name <- paste0("model_FS_", report)
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
  Crop = character(),
  PollDependency = character(),
  AgrIntensity = character(),
  Sites = character(),
  Design = character(),
  Habitat = character(),
  DistanceMeasure = character(),
  MaxDistance = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each report to extract coefficients and store them
for (report in unique_report) {
  # Construct model name
  model_name <- paste0("model_FS_", report)
  dataset_name <- paste0("FS_", report)
  
  # Get the model and dataset from the global environment
  model <- get(model_name)
  data <- get(dataset_name)
  
  # Extract coefficients and their statistics from the model
  coef_summary <- summary(model)$coefficients$cond # specify the conditional model here, i.e. regression part
  
  slope <- coef_summary["log(distance_m + 1)", "Estimate"]
  std_error <- coef_summary["log(distance_m + 1)", "Std. Error"]
  p_value <- coef_summary["log(distance_m + 1)", "Pr(>|z|)"]
  
  # Extract additional information from the dataset
  authors <- unique(data$authors)
  crop <- unique(data$crop)
  agr_intensity <- unique(data$agr_intensity)
  p_dependency <- unique(data$p_dependency)
  sites <- unique(data$sites)
  design <- unique(data$study_design)
  habitat <- unique(data$habitat)
  distance_measure <- unique(data$distance_measure)
  max_distance <- unique(data$max_distance)
  
  # Append the results to the data frame
  results <- rbind(results, data.frame(
    Authors = authors,
    Slope = slope,
    StdError = std_error,
    PValue = p_value,
    Crop = crop,
    PollDependency = p_dependency,
    AgrIntensity = agr_intensity,
    Sites = sites,
    Design = design,
    Habitat = habitat,
    DistanceMeasure = distance_measure,
    MaxDistance = max_distance
  ))
}

# Print the overview of the extracted results
print(results)

# Save the results to a CSV file --> For Stage 2 of the meta-analysis
write.csv(results, "outputs/fruitset/Fruitset-distance stage 1 results.csv", row.names = FALSE)

### Plot model fit ####
# Create folder to save the plots
dir.create(here("outputs", "fruitset", "model fits"), recursive = TRUE, showWarnings = FALSE)

# Loop through each report to generate plots
for (report in unique_report) {
  
  # Construct dataset and model names
  dataset_name <- paste0("FS_", report)
  model_name <- paste0("model_FS_", report)
  
  # Get dataset and model
  dataset <- get(dataset_name)
  model <- get(model_name)
  if (unique(dataset$study_design) == "single distance per site") {model_type <- "(GLM)"
  } else {model_type <- "(GLMM)"
  }
  
  # Create a smooth sequence of distances for prediction
  smooth_distance <- data.frame(distance_m = seq(min(dataset$distance_m), max(dataset$distance_m), length.out = 100))
  
  # Get predictions with ggeffects (population-level, random effects averaged)
  preds <- ggeffects::ggemmeans(model, terms = list(distance_m = smooth_distance$distance_m), 
                                condition = c(repeat_measures = mean(dataset$repeat_measures, na.rm = TRUE)))
  
  # Plot
  p <- ggplot() +
    geom_ribbon(data = preds, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "grey70", alpha = 0.4) +
    geom_line(data = preds, aes(x = x, y = predicted), colour = "blue", linewidth = 1) +
    geom_point(data = dataset, aes(x = distance_m, y = fruitset), colour = "black", alpha = 0.6, size = 3) +
    labs(x = "Distance to natural habitat (m)", y = "Fruit set", title = paste(report, "et al.", model_type)) +
    theme_minimal(base_size = 12)
    
  print(p)
    
  ggsave(here("outputs", "fruitset", "model fits", paste0(report, "_glmmTMB_fit.png")), plot = p, width = 8, height = 6, dpi = 300)
  }
  
########################## Meta-Analysis ###########################

### Conduct the meta-analysis using metafor
# Load the CSV file into a dataframe
fruitset_es <- results
fruitset_es$Variance <- fruitset_es$StdError^2

### Fit a random-effects model using the calculated effect sizes ###
res <- rma(yi = Slope, vi = StdError^2, data = fruitset_es)
res # report tau^2, tau, I^2 and H^2 and test for heterogeinity (Q test)
predict(res, digits=3)
confint(res) # Look at confidence intervals (useful to report this)

# Create the forest plot
tiff(filename = here("outputs", "fruitset", "Fruit set forest plot.tiff"), width = 9, height = 7, units = "in", res = 600)  # Open TIFF device
forest(res,
       slab = fruitset_es$Authors,                # Labels for the studies
       xlab = "Slope",                       # Label for the x-axis
       xlim = c(-2, 2),                            # Customize x-axis limits
       refline = 0,                                # Add reference line at 0
       header = "a) Fruit set",         # Header for the plot
       annotate = TRUE,                            # Add report annotations
       ilab.xpos = -0.015,                         # Adjust position of report effect size labels
       cex = 0.8)                                  # Manage overall font size
dev.off()  # Close device to save the file

##### Plot decay surve
# Define distances (log scale for consistency with Ricketts et al.)
# Get the minimum and maximum distances from the dataset
min_distance <- min(df$distance_m, na.rm = TRUE)
max_distance <- max(df$distance_m, na.rm = TRUE)

# Print the min and max distance values
print(paste("Min distance:", min_distance))
print(paste("Max distance:", max_distance))

# Create a sequence of distances for plotting (e.g., 100 points between min and max)
distance_seq <- seq(min_distance, max_distance, length.out = 100)

# Define mean fruitset at 0m
str(df)
fruitset_at_0m <- df[df$distance_m == 0, ] # Filter the dataset for only rows where distance = 0m
mean_fruitset_0m <- mean(fruitset_at_0m$fruitset, na.rm = TRUE)
print(mean_fruitset_0m) # Print the result
f0 <- mean_fruitset_0m  # mean from dataset = 0.42

# Define slopes from meta-analysis
mean_slope <- res$beta[1]
lower_CI <- res$ci.lb      
upper_CI <- res$ci.ub

# Create a dataframe for the mean effect and confidence intervals
df_fruitset <- data.frame(
  distance_m = distance_seq,
  fruitset_mean = f0 * exp(mean_slope * log(distance_seq + 1)),  # Mean effect
  fruitset_upper = f0 * exp(upper_CI * log(distance_seq + 1)),  # Upper CI
  fruitset_lower = f0 * exp(lower_CI * log(distance_seq + 1))  # Lower CI
)

# Plot using ggplot2
decaycurve <- ggplot(df_fruitset, aes(x = distance_seq)) +
  geom_ribbon(aes(ymin = fruitset_lower, ymax = fruitset_upper), fill = "lightgrey", alpha = 0.5) +  # Shaded CI region
  geom_line(aes(y = fruitset_mean), color = "blue", linewidth = 1) +  # Mean effect
  geom_line(aes(y = fruitset_upper), linetype = "dashed", color = "black") +  # Upper CI
  geom_line(aes(y = fruitset_lower), linetype = "dashed", color = "black") +  # Lower CI
  scale_x_continuous(breaks = seq(0, 5000, by = 1000), limits = c(0, 5000), expand = c(0, 0)) +  # Force x-axis to start at 0
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1), expand = c(0, 0)) +  # Force y-axis to start at 0
  labs(
    x = "Distance to nearest natural habitat (m)",
    y = "Predicted pollinator fruit set (proportion)",
    title = "b) Predicted decay curve") +
  theme_classic(base_size = 14) +  # Standard font size
  theme(
    plot.title = element_text(size = 9, face = "bold"),  # Match forest plot title
    axis.title = element_text(size = 9),  # Match axis labels
    axis.text = element_text(size = 9),  # Standardize tick labels
    plot.margin = margin(10, 20, 10, 10)  # Adjust spacing
  )
ggsave(filename = here("outputs", "fruitset", "Fruitset decay curve.png"), plot = decaycurve, width = 6, height = 5, units = "in", dpi = 300)

#### Calculate % decline at 1km from natural habitat
# Define the function for predicted fruitset
predict_fruitset <- function(distance_m, slope) {
  f0 * exp(slope * log(distance_m + 1))
}

# Calculate predicted fruitsets at 0m and 1000m
predicted_fruitset_0m <- predict_fruitset(0, mean_slope)
predicted_fruitset_1000m <- predict_fruitset(1000, mean_slope)

# Compute percentage decline
percentage_decline <- (1 - (predicted_fruitset_1000m / predicted_fruitset_0m)) * 100

# Print result
print(paste("Predicted % decline at 1km:", round(percentage_decline, 2), "%"))

########################## Additional analyses ###########################

### Moderator analysis ###
# Summarise each moderators distribution
moderators <- c("DistanceMeasure" , "PollDependency", "AgrIntensity", "Habitat") # List of moderator variables

# Use lapply to get counts for each moderator
moderator_summaries <- lapply(moderators, function(moderator) {
  fruitset_es %>%
    count(!!sym(moderator), name = "Study_Count") %>%
    arrange(desc(Study_Count))
})

names(moderator_summaries) <- moderators # Name the list elements
moderator_summaries # Print all summaries

## Meta-analysis with individual moderators
# Moderator analysis for agricultural intensity
res.modintensity <- rma(Slope, Variance, mods = ~ 0 + AgrIntensity, data=fruitset_es)
res.modintensity

# Moderator analysis for habitat type
res.modhabitat <- rma(Slope, Variance, mods = ~ 0 + Habitat, data=fruitset_es)
res.modhabitat

# Moderator analysis for pollinator dependency
res.modpdep <- rma(Slope, Variance, mods = ~ 0 + PollDependency, data=fruitset_es)
res.modpdep

########################## Sensitivity analyses ###########################

# Test for small study bias incl publication bias
metafor::funnel(res) #should be symmetric
regtest(res) # Egger's regression test. If this is not statistically significant, no evidence for small study bias
ranktest(res)

# If small study bias: account for publication bias by using weight function
# wf <- weightr::weightfunct(fruitset_es$Slope, fruitset_es$Variance, table = TRUE)
# wf

### Influential studies ###
inf <- influence(res) # Check if any individual studies were very influential
print(inf) # influential studies have an * next to them
plot(inf) # red dots for influential studies 

# Re-run meta-analysis excluding each report one by one
leave1out <- leave1out(res, digits = 3)
leave1out
write.csv(leave1out, "outputs/fruitset/Leave1out_fruitset.csv", row.names = FALSE)

# QQ-normal plots
qqnorm(res, main = "Random-Effects Model")

# Sensitivity analyses --> methods
res.moddistance <- rma(Slope, Variance, mods = ~ 0 + DistanceMeasure, data=fruitset_es)
res.moddistance

# Add categories for the maximum distance scales (small, medium, large) for sensitivty analysis
fruitset_es <- fruitset_es %>%
  mutate(DistanceCategory = case_when(
    MaxDistance < 750 ~ "small",
    MaxDistance >= 750 & MaxDistance <= 3000 ~ "medium",
    MaxDistance > 3000 ~ "large"
  ))

# Convert DistanceCategory to a factor
fruitset_es$DistanceCategory <- factor(fruitset_es$DistanceCategory, levels = c("small", "medium", "large"))

res.modmaxdistance <- rma(Slope, Variance, mods = ~ 0 + DistanceCategory, data=fruitset_es)
res.modmaxdistance

