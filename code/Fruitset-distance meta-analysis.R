### Two-stage meta-analysis on effects of distance from natural habitat on fruit set

# Load the libraries
library(ggplot2)
library(dplyr)
library(MASS)
library(metafor)
library(here)

# Load data
here() 
df <- read.csv("raw_data/Fulldataset_fruitset.csv", header = TRUE, stringsAsFactors = FALSE)

# Remove rows with NA in 'fruitset' or 'distance' columns
df <- df[!is.na(df$fruitset) & !is.na(df$distance), ]
nrow(df)
colnames(df)

# Get all individual study names
unique_study <- unique(df$study)
unique_study
length(unique_study) #number of studies

# Count unique locations per author and add total sum
author_location_count <- df %>%
  group_by(authors) %>%
  summarise(unique_locations = n_distinct(location)) %>%
  bind_rows(tibble(authors = "Total", unique_locations = sum(.$unique_locations)))
print(author_location_count)

range(df$distance, na.rm = TRUE)

range(df$distance, na.rm = TRUE)

# Calculate the range of distance values for each study
distance_range_by_study <- df %>%
  group_by(study) %>%
  summarise(min_distance = min(distance, na.rm = TRUE),
            max_distance = max(distance, na.rm = TRUE),
            n_observations = n()) # Optional: count number of observations per study

# Print the results
print(distance_range_by_study)

# Calculate summary statistics for max_distance
summary_stats <- distance_range_by_study %>%
  summarise(mean_max_distance = mean(max_distance, na.rm = TRUE),
            median_max_distance = median(max_distance, na.rm = TRUE),
            min_max_distance = min(max_distance, na.rm = TRUE),
            max_max_distance = max(max_distance, na.rm = TRUE))

# Print summary statistics
print(summary_stats)

# Merge max_distance into df
df <- df %>%
  left_join(distance_range_by_study %>% dplyr::select(study, max_distance), by = "study")

# View as a sorted table for better readability
distance_range_by_study %>% arrange(min_distance, max_distance)
write.csv(distance_range_by_study,  "outputs/fruitset/Distance range by study.csv", row.names = FALSE)

########################## Estimate Effect Sizes ###########################

## Estimate effect sizes for each individual study
# Loop through each study and create a separate data frame for each
for (study in unique_study) {
  
  # Filter the df for the current study
  study_data <- df %>% filter(study == !!study)
  
  # Create a variable name dynamically based on the study name
  dataset_name <- paste0("FS_", gsub(" ", "_", study))
  hist(study_data$fruitset)
  
  # Assign the filtered data frame to a new variable in the global environment
  assign(dataset_name, study_data, envir = .GlobalEnv)
}

## Fit GLM Models ###
# Loop through each study to fit the GLM and store the results
for (study in unique_study) {
  # Construct the name of the dataset
  dataset_name <- paste0("FS_", study)
  
  # Access the dataset using get()
  data <- get(dataset_name)
  
  # Fit the GLM
  model <- glm(fruitset ~ log(distance + 1), family = binomial(link = "logit"), data = data)
  
  # Dynamically create the model name in the global environment
  model_name <- paste0("GLM_FS_", study)
  assign(model_name, model, envir = .GlobalEnv)
  
  # Optionally, print the summary of the model
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
  DistanceMeasure = character(),
  MaxDistance = numeric(),
  RoB = character(),
  stringsAsFactors = FALSE
)

# Loop through each study to extract coefficients and store them
for (study in unique_study) {
  # Construct model name
  model_name <- paste0("GLM_FS_", study)
  dataset_name <- paste0("FS_", study)
  
  # Get the model and dataset from the global environment
  model <- get(model_name)
  data <- get(dataset_name)
  
  # Extract coefficients and their statistics from the model
  coef_summary <- summary(model)$coefficients
  
  slope <- coef_summary["log(distance + 1)", "Estimate"]
  std_error <- coef_summary["log(distance + 1)", "Std. Error"]
  p_value <- coef_summary["log(distance + 1)", "Pr(>|z|)"]
  
  # Extract additional information from the dataset
  authors <- unique(data$authors)
  crop <- unique(data$crop)
  agr_intensity <- unique(data$agr_intensity)
  p_dependency <- unique(data$p_dependency)
  sites <- unique(data$sites)
  distance_measure <- unique(data$distance_measure)
  max_distance <- unique(data$max_distance)
  #RoB <- unique(data$RoB)
  
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
    DistanceMeasure = distance_measure,
    MaxDistance = max_distance
    #RoB = RoB
  ))
}

# Print the overview of the extracted results
print(results)

# Save the results to a CSV file --> For Stage 2 of the meta-analysis
write.csv(results, "outputs/fruitset/Fruitset_GLM_Model_Results.csv", row.names = FALSE)

### Plot model fit ####
# Create folder to save the plots
dir.create(here("outputs", "fruitset", "GLM fits"), recursive = TRUE, showWarnings = FALSE)

# Loop through each study to generate plots
for (study in unique_study) {
  # Construct dataset and model names
  dataset_name <- paste0("FS_", study)
  model_name <- paste0("GLM_FS_", study)
  
  # Get dataset and model
  dataset <- get(dataset_name)
  model <- get(model_name)
  
  # Create a smooth sequence of distances for prediction
  smooth_distance <- data.frame(distance = seq(min(dataset$distance), max(dataset$distance), length.out = 100))
  
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
  p <- ggplot(dataset, aes(x = distance, y = fruitset)) +  # Base plot with actual data
    geom_ribbon(data = smooth_distance, aes(x = distance, ymin = lower_CI, ymax = upper_CI), 
                inherit.aes = FALSE, fill = "grey70", alpha = 0.4) +  # Confidence interval shading
    geom_line(data = smooth_distance, aes(x = distance, y = predicted), 
              inherit.aes = FALSE, color = "blue", linewidth = 1) +  # Smooth fitted line
    geom_point(size = 3, alpha = 0.6, colour = "black") + # Raw data points
    labs(x = "Distance in m", y = "Fruit set", title = paste(study, "et al.")) +
    theme_minimal(base_size = 12)  # Clean theme
  
  # Print the plots
  print(p)
  
  # Save the plots
  ggsave(filename = here("outputs", "fruitset", "GLM fits", paste0("Model_Fit_", study, ".png")), 
         plot = p, width = 8, height = 6, dpi = 300)
}

########################## Meta-Analysis ###########################

### Conduct the meta-analysis using metafor
# Load the CSV file into a dataframe
fruitset_es <- read.csv("outputs/fruitset/Fruitset_GLM_Model_Results.csv", header=TRUE, stringsAsFactors = FALSE)
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
       xlim = c(-15, 11),                            # Customize x-axis limits
       refline = 0,                                # Add reference line at 0
       header = "a) Fruit set",         # Header for the plot
       annotate = TRUE,                            # Add study annotations
       ilab.xpos = -0.015,                         # Adjust position of study effect size labels
       cex = 0.8)                                  # Manage overall font size
dev.off()  # Close device to save the file

##### Plot decay surve
# Define distances (log scale for consistency with Ricketts et al.)
# Get the minimum and maximum distances from the dataset
min_distance <- min(df$distance, na.rm = TRUE)
max_distance <- max(df$distance, na.rm = TRUE)

# Print the min and max distance values
print(paste("Min distance:", min_distance))
print(paste("Max distance:", max_distance))

# Create a sequence of distances for plotting (e.g., 100 points between min and max)
distance_seq <- seq(min_distance, max_distance, length.out = 100)

# Define mean fruitset at 0m
str(df)
fruitset_at_0m <- df[df$distance == 0, ] # Filter the dataset for only rows where distance = 0m
mean_fruitset_0m <- mean(fruitset_at_0m$fruitset, na.rm = TRUE)
print(mean_fruitset_0m) # Print the result
f0 <- mean_fruitset_0m  # mean from dataset = 0.42

# Define slopes from meta-analysis
mean_slope <- res$beta[1]
lower_CI <- res$ci.lb      
upper_CI <- res$ci.ub

# Create a dataframe for the mean effect and confidence intervals
df_fruitset <- data.frame(
  distance = distance_seq,
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
predict_fruitset <- function(distance, slope) {
  f0 * exp(slope * log(distance + 1))
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
moderators <- c("DistanceMeasure" , "PollDependency", "AgrIntensity") # List of moderator variables

# Use lapply to get counts for each moderator
moderator_summaries <- lapply(moderators, function(moderator) {
  fruitset_es %>%
    count(!!sym(moderator), name = "Study_Count") %>%
    arrange(desc(Study_Count))
})

names(moderator_summaries) <- moderators # Name the list elements
moderator_summaries # Print all summaries

## Meta-analysis with individual moderators
res.modintensity <- rma(Slope, Variance, mods = ~ 0 + AgrIntensity, data=fruitset_es)
res.modintensity

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

# Re-run meta-analysis excluding each study one by one
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

### Subgroup analysis excl studies with high RoB ###
# Filter for medium RoB studies in the fruitset dataset
# results_mediumRoB <- subset(fruitset_es, RoB == "medium")

# Calculate variance for the meta-analysis
# results_mediumRoB$Variance <- results_mediumRoB$StdError^2

# Fit random-effects model
# res_fruitset_mediumRoB <- rma(yi = Slope, vi = Variance, data = results_mediumRoB)
# print(res_fruitset_mediumRoB)

# Create a forest plot
# forest(res_fruitset_mediumRoB,
#       slab = results_mediumRoB$Authors,                # Labels for the studies
#       xlab = "Slope",                                  # Label for the x-axis
#       xlim = c(-4, 3),                                 # Customise x-axis limits
#       refline = 0,                                     # Add reference line at 0
#       header = "Fruit set subgroup; medium Risk of Bias", # Header for the plot
#       annotate = TRUE,                                 # Add study annotations
#       ilab.xpos = -0.015,                              # Adjust position of study effect size labels
#       cex = 0.8                                        # Manage overall font size
#       )