library(brms)
library(dplyr)
library(bayesplot)
library(ggplot2)
library(loo)
library(posterior)
library(tidyr)
library(broom)
library(patchwork)

# Load in raw data
data <- read.csv("ewpw_chick_measurement_data.csv") 

# Filter data to only include mass and tarsus measurements
mass_tarsus_data <- data %>%
  dplyr::select(NestId, BandId, age_days, mass, tarsus)

# Filter data to include all other measurements and exclude the first two nests because the other measurements were not taken at these nests
other_data <- data %>%
  filter(!NestId %in% c('EPPLNest1', 'EPPLNest2')) %>%
  dplyr::select(NestId, BandId, age_days, bill_length_from_nares, posterior_bill_width, anterior_bill_width, 
                mouth_width, humeral_length, tail_length, wing_length)

################################################################################

# Use wide diffuse priors for the Gompertz model
priors <- c(
  prior(normal(0, 100), nlpar = "a"),
  prior(normal(0, 100), nlpar = "b"),
  prior(normal(0, 100), nlpar = "c"))

# Define the Gompertz equation for each morphometric measurement and add the priors and random effect for individual
formula_mass <- bf(mass ~ a * exp(-b * exp(-c * age_days)), 
                   a ~ 1 + (1|BandId), 
                   b ~ 1 + (1|BandId), 
                   c ~ 1 + (1|BandId), 
                   nl = TRUE)

formula_tarsus <- bf(tarsus ~ a * exp(-b * exp(-c * age_days)), 
                     a ~ 1 + (1|BandId), 
                     b ~ 1 + (1|BandId), 
                     c ~ 1 + (1|BandId), 
                     nl = TRUE)

formula_humeral_length <- bf(humeral_length ~ a * exp(-b * exp(-c * age_days)), 
                             a ~ 1 + (1|BandId), 
                             b ~ 1 + (1|BandId), 
                             c ~ 1 + (1|BandId), 
                             nl = TRUE)

formula_wing_length <- bf(wing_length ~ a * exp(-b * exp(-c * age_days)), 
                          a ~ 1 + (1|BandId), 
                          b ~ 1 + (1|BandId), 
                          c ~ 1 + (1|BandId), 
                          nl = TRUE)

formula_bill_length_from_nares <- bf(bill_length_from_nares ~ a * exp(-b * exp(-c * age_days)), 
                                     a ~ 1 + (1|BandId), 
                                     b ~ 1 + (1|BandId), 
                                     c ~ 1 + (1|BandId), 
                                     nl = TRUE)

formula_posterior_bill_width <- bf(posterior_bill_width ~ a * exp(-b * exp(-c * age_days)), 
                                   a ~ 1 + (1|BandId), 
                                   b ~ 1 + (1|BandId), 
                                   c ~ 1 + (1|BandId), 
                                   nl = TRUE)

formula_anterior_bill_width <- bf(anterior_bill_width ~ a * exp(-b * exp(-c * age_days)), 
                                  a ~ 1 + (1|BandId), 
                                  b ~ 1 + (1|BandId), 
                                  c ~ 1 + (1|BandId), 
                                  nl = TRUE)

formula_mouth_width <- bf(mouth_width ~ a * exp(-b * exp(-c * age_days)), # Mouth width would not fit when including the random effect so we removed it, this was noted in the manuscript
                          a ~ 1 , 
                          b ~ 1 , 
                          c ~ 1 , 
                          nl = TRUE)

formula_tail_length <- bf(tail_length ~ a * exp(-b * exp(-c * age_days)), 
                          a ~ 1 + (1|BandId), 
                          b ~ 1 + (1|BandId), 
                          c ~ 1 + (1|BandId), 
                          nl = TRUE)

################################################################################

# Fit the simulations for each of the morphometric measurements
fit_mass <- brm(
  formula = formula_mass,
  data = mass_tarsus_data,
  family = gaussian(),
  prior = priors,
  seed = 12345,
  chains = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 14,
  save_pars = save_pars(all = TRUE)
)

fit_tarsus <- brm(
  formula = formula_tarsus,
  data = mass_tarsus_data,
  family = gaussian(),
  prior = priors,
  seed = 12345,
  chains = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 14,
  save_pars = save_pars(all = TRUE)
)

fit_humeral_length <- brm(
  formula = formula_humeral_length,
  data = other_data,
  family = gaussian(),
  prior = priors,
  seed = 12345,
  chains = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 14,
  save_pars = save_pars(all = TRUE)
)      

fit_wing_length <- brm(
  formula = formula_wing_length,
  data = other_data,
  family = gaussian(),
  prior = priors,
  seed = 12345,
  chains = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 14,
  save_pars = save_pars(all = TRUE)
)      

fit_bill_length_from_nares <- brm(
  formula = formula_bill_length_from_nares,
  data = other_data,
  family = gaussian(),
  prior = priors,
  seed = 12345,
  chains = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 14,
  save_pars = save_pars(all = TRUE)
)

fit_posterior_bill_width <- brm(
  formula = formula_posterior_bill_width,
  data = other_data,
  family = gaussian(),
  prior = priors,
  seed = 12345,
  chains = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 14,
  save_pars = save_pars(all = TRUE)
)                

fit_anterior_bill_width <- brm(
  formula = formula_anterior_bill_width,
  data = other_data,
  family = gaussian(),
  prior = priors,
  seed = 12345,
  chains = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 14,
  save_pars = save_pars(all = TRUE)
)    

fit_mouth_width <- brm(
  formula = formula_mouth_width,
  data = other_data,
  family = gaussian(),
  prior = priors,
  seed = 12345,
  chains = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 14,
  save_pars = save_pars(all = TRUE)
)                

fit_tail_length <- brm(
  formula = formula_tail_length,
  data = other_data,
  family = gaussian(),
  prior = priors,
  seed = 12345,
  chains = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  cores = 14,
  save_pars = save_pars(all = TRUE)
)      

# Get summary stats outputs for each model
summary(fit_mass)
summary(fit_tarsus)
summary(fit_humeral_length)
summary(fit_wing_length)
summary(fit_bill_length_from_nares)
summary(fit_posterior_bill_width)
summary(fit_anterior_bill_width)
summary(fit_mouth_width)
summary(fit_tail_length)

# Save models 
# saveRDS(fit_mass, "fit_mass.rds")
# saveRDS(fit_tarsus, "fit_tarsus.rds")
# saveRDS(fit_humeral_length, "fit_humeral_length.rds")
# saveRDS(fit_wing_length, "fit_wing_length.rds")
# saveRDS(fit_bill_length_from_nares, "fit_bill_length_from_nares.rds")
# saveRDS(fit_posterior_bill_width, "fit_posterior_bill_width.rds")
# saveRDS(fit_anterior_bill_width, "fit_anterior_bill_width.rds")
# saveRDS(fit_mouth_width, "fit_mouth_width.rds")
# saveRDS(fit_tail_length, "fit_tail_length.rds")


# Run MCMC trace plots
mcmc_trace(fit_mass)
mcmc_trace(fit_tarsus)
mcmc_trace(fit_humeral_length)
mcmc_trace(fit_wing_length)
mcmc_trace(fit_bill_length_from_nares)
mcmc_trace(fit_posterior_bill_width)
mcmc_trace(fit_anterior_bill_width)
mcmc_trace(fit_mouth_width)
mcmc_trace(fit_tail_length)

# Run Posterior Predictive Checks
pp_check(fit_mass)
pp_check(fit_tarsus)
pp_check(fit_humeral_length)
pp_check(fit_wing_length)
pp_check(fit_bill_length_from_nares)
pp_check(fit_posterior_bill_width)
pp_check(fit_anterior_bill_width)
pp_check(fit_mouth_width)
pp_check(fit_tail_length)

# Compute LOO-adjusted R-squared
loo_R2(fit_mass)
loo_R2(fit_tarsus)
loo_R2(fit_humeral_length)
loo_R2(fit_wing_length)
loo_R2(fit_bill_length_from_nares)
loo_R2(fit_posterior_bill_width)
loo_R2(fit_anterior_bill_width)
loo_R2(fit_mouth_width)
loo_R2(fit_tail_length)

# Identify the most informative models 

################################################################################

# Load saved models if needed
fit_mass <- readRDS("fit_mass.rds")
fit_tarsus <- readRDS("fit_tarsus.rds")
fit_humeral_length <- readRDS("fit_humeral_length.rds")
fit_wing_length <- readRDS("fit_wing_length.rds")

# Create growth curves

# Define a sequence of ages for prediction
age_seq <- seq(from = min(other_data$age_days), 
               to = max(other_data$age_days), 
               length.out = 100)

# Extract a unique BandId for use in plotting
unique_band_id <- unique(fit_mass$data$BandId)[1]  

# Function to generate new data and run predictions
predict_growth <- function(model, variable_name) {
  # Extract a valid BandId from the model data
  unique_band_id <- unique(model$data$BandId)[1]  
  
  # Create a new data frame for predictions
  new_data <- data.frame(
    age_days = seq(from = min(model$data$age_days), 
                   to = max(model$data$age_days), 
                   length.out = 100),
    BandId = unique_band_id 
  )
  
  # Generate posterior predictions
  posterior_samples <- posterior_predict(model, newdata = new_data, allow_new_levels = TRUE)
  
  # Compute median and credible intervals
  pred_df <- tibble(
    age_days = new_data$age_days,
    mean = apply(posterior_samples, 2, mean),
    lower = apply(posterior_samples, 2, quantile, probs = 0.025),
    upper = apply(posterior_samples, 2, quantile, probs = 0.975)
  )
  
  return(pred_df)
}

# Run predictions for all four models
pred_mass <- predict_growth(fit_mass, "Mass")
pred_tarsus <- predict_growth(fit_tarsus, "Tarsus")
pred_humeral <- predict_growth(fit_humeral_length, "Humeral Length")
pred_wing <- predict_growth(fit_wing_length, "Wing Length")

# Combine all predictions into a single data frame
pred_mass$Measurement <- "Mass"
pred_tarsus$Measurement <- "Tarsus"
pred_humeral$Measurement <- "Humerus"
pred_wing$Measurement <- "Wing"

pred_data_all <- bind_rows(pred_mass, pred_tarsus, pred_humeral, pred_wing)

# Plot the growth curves
p <- ggplot(pred_data_all %>% filter(!is.na(mean)),  
            aes(x = age_days, y = mean, group = Measurement)) +
  geom_line(aes(linetype = Measurement), size = 1, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.4) +
  
  # Labels for each measurement
  annotate("text", x = 15, y = 18, label = "Tarsus",  
           hjust = -0.1, size = 5.5, family = "serif", color = "black") +
  annotate("text", x = 15, y = 85, label = "Wing",  
           hjust = -0.18, size = 5.5, family = "serif", color = "black") +
  annotate("text", x = 15, y = 34, label = "Humerus",  
           hjust = -0.1, size = 5.5, family = "serif", color = "black") +
  annotate("text", x = 15, y = 41, label = "Mass",  
           hjust = -0.18, size = 5.5, family = "serif", color = "black") +
  annotate("text", x = 0, y = 95, label = "a",  
           hjust = -0.5, size = 8, family = "serif", color = "black") +
  
  # Labels and theme settings
  labs(x = "Days", y = "Measurement value (g or mm)") +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    text = element_text(family = "serif", size = 16),  
    axis.title.x = element_text(size = 16, margin = margin(t = 12), color = "black"),  
    axis.title.y = element_text(size = 16, margin = margin(r = 12), color = "black"),  
    axis.text.x = element_text(size = 16, color = "black"),  
    axis.text.y = element_text(size = 16, color = "black"),  
    plot.margin = unit(c(1, 3, 1, 1), "cm"),  
    legend.position = "none"
  ) +
  coord_cartesian(clip = "off")

# Print the plot
print(p)

################################################################################

# Create growth rate curves

# Gompertz Growth Rate Function
gompertz_derivative <- function(t, a, b, c) {
  a * b * c * exp(-b * exp(-c * t)) * exp(-c * t)
}

# Function to compute growth rates from posterior samples
compute_growth_rate_gompertz <- function(model, variable_name) {
  # Extract posterior samples for parameters a, b, and c
  post_samples <- as_draws_df(model)
  
  # Ensure correct parameter names based on model structure
  a_samples <- post_samples[[paste0("b_a_Intercept")]]
  b_samples <- post_samples[[paste0("b_b_Intercept")]]
  c_samples <- post_samples[[paste0("b_c_Intercept")]]
  
  # Define time sequence (Start at Day 1 instead of 0)
  age_seq <- seq(from = 1, 
                 to = max(model$data$age_days), 
                 length.out = 100)
  
  # Compute growth rates for each posterior draw
  growth_rates <- sapply(1:length(a_samples), function(i) {
    gompertz_derivative(age_seq, a_samples[i], b_samples[i], c_samples[i])
  })
  
  # Compute median and credible intervals
  growth_rate_df <- tibble(
    age_days = age_seq,
    mean = apply(growth_rates, 1, median),
    lower = apply(growth_rates, 1, quantile, probs = 0.025),
    upper = apply(growth_rates, 1, quantile, probs = 0.975),
    Measurement = variable_name
  )
  
  return(growth_rate_df)
}

# Calculate growth rates for each measurement
growth_mass <- compute_growth_rate_gompertz(fit_mass, "Mass")
growth_tarsus <- compute_growth_rate_gompertz(fit_tarsus, "Tarsus")
growth_humeral <- compute_growth_rate_gompertz(fit_humeral_length, "Humerus")
growth_wing <- compute_growth_rate_gompertz(fit_wing_length, "Wing")

# Combine all growth rate data into one data frame
growth_data_all <- bind_rows(growth_mass, growth_tarsus, growth_humeral, growth_wing)

# Ensure x-axis starts at 0 (to match growth curve plot) but remove growth rate values at day 0 (growth rates can't start on day 0)
growth_data_all <- growth_data_all %>% 
  mutate(mean = ifelse(age_days == 0, NA, mean),
         lower = ifelse(age_days == 0, NA, lower),
         upper = ifelse(age_days == 0, NA, upper))

# Plot the growth rates
p_growth_rates <- ggplot(growth_data_all, aes(x = age_days, y = mean, group = Measurement)) +
  geom_line(aes(linetype = Measurement), size = 1, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.4) +
  
  # Add labels for each measurement (adjust y-values as needed)
  annotate("text", x = 15, y = 0.13, label = "Tarsus",  
           hjust = -0.15, size = 5.5, family = "serif", color = "black") +
  annotate("text", x = 15, y = 4.2, label = "Wing",  
           hjust = -0.1, size = 5.5, family = "serif", color = "black") +
  annotate("text", x = 15, y = 0.55, label = "Humerus",  
           hjust = -0.1, size = 5.5, family = "serif", color = "black") +
  annotate("text", x = 15, y = 1.1, label = "Mass",  
           hjust = -0.2, size = 5.5, family = "serif", color = "black") +
  annotate("text", x = 0, y = 6.4, label = "b",  
           hjust = -0.5, size = 8, family = "serif", color = "black") +  
  
  # Labels and theme settings
  labs(x = "Days", y = "Growth rate (g or mm per day)") +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    text = element_text(family = "serif", size = 16),  
    axis.title.x = element_text(size = 16, margin = margin(t = 12), color = "black"),  
    axis.title.y = element_text(size = 16, margin = margin(r = 12), color = "black"),  
    axis.text.x = element_text(size = 16, color = "black"),  
    axis.text.y = element_text(size = 16, color = "black"),  
    plot.margin = unit(c(1, 5, 1, 1), "cm"),  
    legend.position = "none"
  ) +
  coord_cartesian(clip = "off", xlim = c(0, max(growth_data_all$age_days)))  # Ensures x-axis starts at 0

# Create a combined plot
combined_plot <- p + p_growth_rates + plot_layout(ncol = 1)

print(combined_plot)

# Save the combined figure
# ggsave(filename = "D:/Figures/growth_curve_and_rates_updated.png", 
       # plot = combined_plot, 
       # width = 8, height = 10, dpi = 300, units = "in")

################################################################################

# Identify max growth rate for each variable with whole number days (1-15)
max_growth_rates_discrete <- growth_data_all %>%
  mutate(age_days = round(age_days)) %>%   # Round age days to nearest whole number
  filter(age_days >= 1, age_days <= 15) %>%  # Keep only whole days 1-15
  group_by(Measurement, age_days) %>%  
  summarise(mean_growth = mean(mean, na.rm = TRUE), .groups = "drop") %>%  # Average if multiple rows exist per day
  group_by(Measurement) %>%
  summarise(
    max_growth_rate = max(mean_growth, na.rm = TRUE),  # Find max growth rate
    day_at_max = age_days[which.max(mean_growth)]      # Find corresponding whole day
  )

# Print results
print(max_growth_rates_discrete)

################################################################################

#Calculate growth rates with methods in Ramstack et al. 1998 for comparison to swallow species

# Function to compute growth rates using linear regression
compute_growth_rate_lm <- function(data, measurement, day_min, day_max) {
  # Filter data for the specified age range
  filtered_data <- data %>%
    filter(age_days >= day_min, age_days <= day_max) 
  
  # Run linear model (Measurement ~ Age)
  growth_model <- lm(as.formula(paste(measurement, "~ age_days")), data = filtered_data)
  
  # Extract growth rate
  summary_df <- tidy(growth_model) %>%
    filter(term == "age_days") %>% 
    select(estimate) %>%
    rename(growth_rate = estimate) %>%
    mutate(Measurement = measurement)
  
  return(summary_df)
}

# Compute growth rates for each measurement
growth_rates_mass <- compute_growth_rate_lm(mass_tarsus_data, "mass", 3, 10)
growth_rates_tarsus <- compute_growth_rate_lm(mass_tarsus_data, "tarsus", 1, 10)
growth_rates_wing <- compute_growth_rate_lm(other_data, "wing_length", 5, max(other_data$age_days))

# Combine results into one table
growth_rates_summary <- bind_rows(growth_rates_mass, growth_rates_tarsus, growth_rates_wing)

print(growth_rates_summary)

################################################################################

# Create the predictive tables for each model

new_data <- data.frame(age_days = 0:15, BandId = factor(rep("new_band", 16)))

# Function to generate predictions for a given model and measurement name
generate_prediction_table <- function(model, variable_name) {
  # Generate posterior predictions including random effects
  predicted_variable <- posterior_epred(model, newdata = new_data, allow_new_levels = TRUE)
  
  # Calculate mean and 95% credible intervals
  variable_summary <- apply(predicted_variable, 2, function(x) {
    c(mean = mean(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
  })
  
  # Convert to dataframe
  variable_summary_df <- as.data.frame(t(variable_summary))
  variable_summary_df$Age_Days <- new_data$age_days
  
  # Rename columns
  names(variable_summary_df) <- c("Predicted_Mean", "Lower_CI", "Upper_CI", "Age_Days")
  
  # Add measurement type
  variable_summary_df$Measurement <- variable_name
  
  # Round numeric values to 1 decimal place
  variable_summary_df <- variable_summary_df %>%
    mutate(across(where(is.numeric), ~ round(.x, 1)))
  
  return(variable_summary_df)
}

# Run the function for all four models
pred_mass <- generate_prediction_table(fit_mass, "Mass")
pred_tarsus <- generate_prediction_table(fit_tarsus, "Tarsus")
pred_wing <- generate_prediction_table(fit_wing_length, "Wing")
pred_humerus <- generate_prediction_table(fit_humeral_length, "Humerus")

# Combine results into one table
predictions_all <- bind_rows(pred_mass, pred_tarsus, pred_wing, pred_humerus)

# Reorder columns
predictions_all <- predictions_all %>% select(Age_Days, Measurement, Predicted_Mean, Lower_CI, Upper_CI)

# Print final table
print(predictions_all)

# Combine all 4 variables into a single unitless value

# Generate predictive data
new_data <- data.frame(age_days = 0:15, BandId = factor(rep("new_band", 16)))

# Generate posterior predictions for each variable
predicted_mass <- posterior_epred(fit_mass, newdata = new_data, allow_new_levels = TRUE)
predicted_tarsus <- posterior_epred(fit_tarsus, newdata = new_data, allow_new_levels = TRUE)
predicted_humerus <- posterior_epred(fit_humeral_length, newdata = new_data, allow_new_levels = TRUE)
predicted_wing <- posterior_epred(fit_wing_length, newdata = new_data, allow_new_levels = TRUE)

# Sum the predictions across the four models and average them
combined_predictions <- (predicted_mass + predicted_tarsus + predicted_humerus + predicted_wing) / 4

# Compute mean and 95% credible intervals for each day
summary_stats <- apply(combined_predictions, 2, function(column) {
  c(
    Mean = mean(column),
    Lower_CI = quantile(column, probs = 0.025),
    Upper_CI = quantile(column, probs = 0.975)
  )
})

# Convert to a data frame
summary_stats_df <- as.data.frame(t(summary_stats))

# Add age_days column
summary_stats_df$Age_Days <- new_data$age_days

# Rename columns 
names(summary_stats_df) <- c("Mean", "Lower_CI", "Upper_CI", "Age_Days")

# Reorder columns correctly
summary_stats_df <- summary_stats_df %>% select(Age_Days, Mean, Lower_CI, Upper_CI)

# Round numeric values to 1 decimal place
summary_stats_df <- summary_stats_df %>% mutate(across(where(is.numeric), ~ round(.x, 1)))

# Print final table
print(summary_stats_df)

################################################################################

# Create growth curves for the remaining variables 

# Load saved models if needed
# fit_bill_length_from_nares <- readRDS("fit_bill_length_from_nares.rds")
# fit_posterior_bill_width <- readRDS("fit_posterior_bill_width.rds")
# fit_anterior_bill_width <- readRDS("fit_anterior_bill_width.rds")
# fit_mouth_width <- readRDS("fit_mouth_width.rds")
# fit_tail_length <- readRDS("fit_tail_length.rds")

# Create a sequence of ages for prediction
age_seq <- seq(from = min(other_data$age_days), 
               to = max(other_data$age_days), 
               length.out = 100)

# Generate a new data frame for predictions
new_data <- data.frame(age_days = age_seq)

# Function to generate predictions for each measurement 
generate_predictions <- function(fit_model, measurement_name) {
  posterior_draws <- posterior_predict(fit_model, newdata = new_data, allow_new_levels = TRUE)
  
  # Compute mean and 95% CI
  pred_summary <- apply(posterior_draws, 2, function(x) {
    c(mean = mean(x), ci_lower = quantile(x, probs = 0.025), ci_upper = quantile(x, probs = 0.975))
  })
  
  pred_data <- as.data.frame(t(pred_summary))
  pred_data$age_days <- age_seq
  pred_data$Measurement <- measurement_name
  
  return(pred_data)
}

# Generate predictions for each measurement 
pred_bill_nares <- generate_predictions(fit_bill_length_from_nares, "Bill Length from Nares")
pred_post_bill <- generate_predictions(fit_posterior_bill_width, "Posterior Bill Width")
pred_ant_bill <- generate_predictions(fit_anterior_bill_width, "Anterior Bill Width")
pred_mouth <- generate_predictions(fit_mouth_width, "Mouth Width")
pred_tail <- generate_predictions(fit_tail_length, "Tail Length")

# Combine all predictions into one data frame
combined_pred_data <- bind_rows(pred_bill_nares, pred_post_bill, pred_ant_bill, pred_mouth, pred_tail)

# Rename columns 
colnames(combined_pred_data) <- c("mean", "ci_lower", "ci_upper", "age_days", "Measurement")

# Reshape other_data to long format 
observed_data_long <- other_data %>%
  pivot_longer(
    cols = c(bill_length_from_nares, posterior_bill_width, anterior_bill_width, mouth_width, tail_length), 
    names_to = "Measurement", values_to = "value"
  ) %>%
  mutate(Measurement = recode(Measurement,
                              bill_length_from_nares = "Bill Length from Nares",
                              posterior_bill_width = "Posterior Bill Width",
                              anterior_bill_width = "Anterior Bill Width",
                              mouth_width = "Mouth Width",
                              tail_length = "Tail Length"))

# Generate the plot 
combined_5_measurements <- ggplot(combined_pred_data, aes(x = age_days, y = mean, shape = Measurement)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.4, fill = "gray") +  
  geom_line(size = 1, colour = "black", linetype = "solid") + 
  geom_point(data = observed_data_long, aes(x = age_days, y = value, shape = Measurement), 
             colour = "black", alpha = 0.7, size = 2.5) +  
  labs(
    x = "Days",
    y = "Millimeters",
    shape = "Measurement" 
  ) +
  facet_wrap(~ Measurement, scales = "free_y") +  
  scale_shape_manual(values = c(
    "Bill Length from Nares" = 16,  # Circle
    "Posterior Bill Width" = 17,    # Triangle
    "Anterior Bill Width" = 15,     # Square
    "Mouth Width" = 18,             # Diamond
    "Tail Length" = 3               # Plus
  )) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    text = element_text(family = "serif", size = 16),  
    axis.title.x = element_text(size = 16, margin = margin(t = 12), color = "black"),  
    axis.title.y = element_text(size = 16, margin = margin(r = 12), color = "black"),  
    axis.text.x = element_text(size = 16, color = "black"),  
    axis.text.y = element_text(size = 16, color = "black"),  
    plot.margin = unit(c(1, 3, 1, 1), "cm"),  
    legend.text = element_text(size = 16),   
    legend.title = element_text(size = 16), 
    legend.position = "right",  
    strip.text = element_blank()  
  )

# Print the updated plot
print(combined_5_measurements)

# Save plot 
# ggsave(
#   filename = "D:/Figures/other_5_measurements_updated.jpg",  
#   plot = combined_5_measurements,        
#   dpi = 300,                             
#   height = 10,                           
#   width = 18,                           
#   units = "in",                          
#   device = "jpeg",                      
#   quality = 100                          
# )
