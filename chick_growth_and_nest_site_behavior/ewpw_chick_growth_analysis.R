# Install packages if needed

# install.packages("brms")
# install.packages("dplyr")
# install.packages("bayesplot")
# install.packages("ggplot2")
# install.packages("loo")
# install.packages("posterior")
# install.packages("tidyr")
# install.packages("patchwork")

# Load in packages 

library(brms)
library(dplyr)
library(bayesplot)
library(ggplot2)
library(loo)
library(posterior)
library(tidyr)
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

######################################################################################################################

# Use wide uninformative priors for the Gompertz model
priors <- c(
  prior(normal(0, 100), nlpar = "a"),
  prior(normal(0, 100), nlpar = "b"),
  prior(normal(0, 100), nlpar = "c"))

# Define the Gompertz equation for each morphometric measurement and add the priors and radon effect for individual
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

######################################################################################################################

# Fit the models for each of the morphometric measurements
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

######################################################################################################################

# MCMC diagnostics 
mcmc_trace(fit_mass) # Replace fit_mass with each of the morphometric measurements
mcmc_pairs(fit_mass) # Replace fit_mass with each of the morphometric measurements

# Posterior predictive checks
pp_check(fit_mass) # Replace fit_mass with each of the morphometric measurements

#Leave-one-out cross validation
loo_result <- loo(fit_mass, moment_match = TRUE, cores = 14) # Replace fit_mass with each of the morphometric measurements

# Check diagnostics after moment matching
print(loo_result)

# Optionally save model outputs
# saveRDS(fit_mass, "fit_mass.rds")
# saveRDS(fit_tarsus, "fit_tarsus.rds")
# saveRDS(fit_humeral_length, "fit_humeral_length.rds")
# saveRDS(fit_wing_length, "fit_wing_length.rds")
# saveRDS(fit_bill_length_from_nares, "fit_bill_length_from_nares.rds")
# saveRDS(fit_posterior_bill_width, "fit_posterior_bill_width.rds")
# saveRDS(summary(fit_anterior_bill_width), "summary(fit_anterior_bill_width).rds")
# saveRDS(fit_mouth_width, "fit_mouth_width.rds")
# saveRDS(fit_tail_length, "fit_tail_length.rds")

######################################################################################################################

# Calculate coefficient of variance for each model
calculate_cv <- function(fit) {
  model_summary <- summary(fit)
  pop_effects <- model_summary$fixed
  estimates <- pop_effects[, "Estimate"]
  std_devs <- pop_effects[, "Est.Error"]
  cv_values <- (std_devs / estimates) * 100
  names(cv_values) <- paste("CV", rownames(pop_effects), sep="_")
  return(cv_values)
}

cv_mass <- calculate_cv(fit_mass)
cv_tarsus <- calculate_cv(fit_tarsus)
cv_fit_humeral_length <- calculate_cv(fit_humeral_length)
cv_fit_wing_length <- calculate_cv(fit_wing_length)
cv_bill_length_from_nares <- calculate_cv(fit_bill_length_from_nares)
cv_fit_posterior_bill_width <- calculate_cv(fit_posterior_bill_width)
cv_fit_anterior_bill_width <- calculate_cv(fit_anterior_bill_width)
cv_fit_tail_length <- calculate_cv(fit_tail_length)
cv_fit_mouth_width <- calculate_cv(fit_mouth_width)

cv_values <- list(
  Mass = calculate_cv(fit_mass),
  Tarsus = calculate_cv(fit_tarsus),
  Humeral_Length = calculate_cv(fit_humeral_length),
  Wing_Length = calculate_cv(fit_wing_length),
  Bill_Length_from_Nares = calculate_cv(fit_bill_length_from_nares),
  Posterior_Bill_Width = calculate_cv(fit_posterior_bill_width),
  Anterior_Bill_Width = calculate_cv(fit_anterior_bill_width),
  Tail_Length = calculate_cv(fit_tail_length),
  Mouth_Width = calculate_cv(fit_mouth_width)
)

# Create a data frame 
cv_data_frame <- do.call(cbind, cv_values)
cv_data_frame <- data.frame(Model = rownames(cv_data_frame), cv_data_frame, row.names = NULL)

# Print the resulting data frame
print(cv_data_frame)

######################################################################################################################

# Generate posterior predictions for the top 4 variables (CV less than 15%): mass, tarsus, humeral length, wing length

# Load saved models if needed
fit_mass <- readRDS("fit_mass.rds")
fit_tarsus <- readRDS("fit_tarsus.rds")
fit_humeral_length <- readRDS("fit_humeral_length.rds")
fit_wing_length <- readRDS("fit_wing_length.rds")


# Define a sequence of ages for prediction
age_seq <- seq(from = min(other_data$age_days), 
               to = max(other_data$age_days), 
               length.out = 100)

# Create a new data frame for the prediction range
new_data <- data.frame(age_days = age_seq)

# Calculate predictions for each model including credible intervals
pred_mass <- posterior_predict(fit_mass, newdata = new_data, re_formula = NA)
pred_tarsus <- posterior_predict(fit_tarsus, newdata = new_data, re_formula = NA)
pred_humeral_length <- posterior_predict(fit_humeral_length, newdata = new_data, re_formula = NA)
pred_wing_length <- posterior_predict(fit_wing_length, newdata = new_data, re_formula = NA)

# Calculate mean and credible intervals for each prediction
summary_data <- function(pred_matrix, age_seq) {
  tibble(
    age_days = age_seq,
    Mean = apply(pred_matrix, 2, mean),
    Lower = apply(pred_matrix, 2, quantile, probs = 0.025),
    Upper = apply(pred_matrix, 2, quantile, probs = 0.975)
  )
}

# Apply the summary function to each prediction matrix
summary_mass <- summary_data(pred_mass, age_seq)
summary_tarsus <- summary_data(pred_tarsus, age_seq)
summary_humeral_length <- summary_data(pred_humeral_length, age_seq)
summary_wing_length <- summary_data(pred_wing_length, age_seq)

# Combine all summaries into one data frame for plotting
pred_data_all <- bind_rows(
  summary_mass %>% mutate(Measurement = "Mass"),
  summary_tarsus %>% mutate(Measurement = "Tarsus"),
  summary_humeral_length %>% mutate(Measurement = "Humeral Length"),
  summary_wing_length %>% mutate(Measurement = "Wing Length")
)

# Define the desired order of the panels
pred_data_all$Measurement <- factor(pred_data_all$Measurement, levels = c("Mass", "Tarsus", "Humeral Length", "Wing Length"))

# Calculate max growth rates for the top 4 variables

# Extract posterior samples for parameters from the fitted models
posterior_samples_mass <- as_draws_df(fit_mass)
posterior_samples_tarsus <- as_draws_df(fit_tarsus)
posterior_samples_humeral_length <- as_draws_df(fit_humeral_length)
posterior_samples_wing_length <- as_draws_df(fit_wing_length)

# Extracting parameter samples for 'a', 'b', and 'c'
posterior_samples_mass <- data.frame(
  a = posterior_samples_mass$b_a_Intercept,
  b = posterior_samples_mass$b_b_Intercept,
  c = posterior_samples_mass$b_c_Intercept
)

posterior_samples_tarsus <- data.frame(
  a = posterior_samples_tarsus$b_a_Intercept,
  b = posterior_samples_tarsus$b_b_Intercept,
  c = posterior_samples_tarsus$b_c_Intercept
)

posterior_samples_humeral_length <- data.frame(
  a = posterior_samples_humeral_length$b_a_Intercept,
  b = posterior_samples_humeral_length$b_b_Intercept,
  c = posterior_samples_humeral_length$b_c_Intercept
)

posterior_samples_wing_length <- data.frame(
  a = posterior_samples_wing_length$b_a_Intercept,
  b = posterior_samples_wing_length$b_b_Intercept,
  c = posterior_samples_wing_length$b_c_Intercept
)


# Define the first derivative of the Gompertz equation
gompertz_derivative <- function(t, a, b, c) {
  a * b * c * exp(-b * exp(-c * t)) * exp(-c * t)
}

# Function to calculate and summarize growth rates
calculate_growth_rate_summary <- function(posterior_samples) {
  # Generate a sequence of time points from 1 to 15
  time_points <- 1:15
  
  # Initialize a list to store the growth rate samples for each time point
  growth_rate_samples <- vector("list", length(time_points))
  
  # Evaluate the derivative for each time point using the posterior samples
  for (i in seq_along(time_points)) {
    t <- time_points[i]
    growth_rate_samples[[i]] <- apply(posterior_samples, 1, function(params) {
      gompertz_derivative(t, params["a"], params["b"], params["c"])
    })
  }
  
  # Convert the list of growth rate samples into a data frame
  growth_rate_df <- do.call(cbind, growth_rate_samples)
  colnames(growth_rate_df) <- paste0("Day_", time_points)
  
  # Summarize the growth rate distribution for each time point
  growth_rate_summary <- data.frame(
    Day = time_points,
    Mean = apply(growth_rate_df, 2, mean),
    Median = apply(growth_rate_df, 2, median),
    Lower_95 = apply(growth_rate_df, 2, quantile, 0.025),
    Upper_95 = apply(growth_rate_df, 2, quantile, 0.975)
  )
  
  return(growth_rate_summary)
}

# Calculate growth rate summaries for each dataset
growth_rate_summary_mass <- calculate_growth_rate_summary(posterior_samples_mass)
growth_rate_summary_tarsus <- calculate_growth_rate_summary(posterior_samples_tarsus)
growth_rate_summary_humeral_length <- calculate_growth_rate_summary(posterior_samples_humeral_length)
growth_rate_summary_wing_length <- calculate_growth_rate_summary(posterior_samples_wing_length)

# Print the summaries
print(growth_rate_summary_mass)
print(growth_rate_summary_tarsus)
print(growth_rate_summary_humeral_length)
print(growth_rate_summary_wing_length)

# Function to find the maximum growth rate and corresponding day
find_max_growth_rate <- function(growth_rate_summary) {
  max_index <- which.max(growth_rate_summary$Mean)
  max_day <- growth_rate_summary$Day[max_index]
  max_growth_rate <- growth_rate_summary$Mean[max_index]
  list(max_day = max_day, max_growth_rate = max_growth_rate)
}

# Find maximum growth rate and corresponding day for each dataset
max_growth_mass <- find_max_growth_rate(growth_rate_summary_mass)
max_growth_tarsus <- find_max_growth_rate(growth_rate_summary_tarsus)
max_growth_humeral_length <- find_max_growth_rate(growth_rate_summary_humeral_length)
max_growth_wing_length <- find_max_growth_rate(growth_rate_summary_wing_length)

# Print the results
print(paste("Mass: Maximum growth rate is", max_growth_mass$max_growth_rate, "on day", max_growth_mass$max_day))
print(paste("Tarsus: Maximum growth rate is", max_growth_tarsus$max_growth_rate, "on day", max_growth_tarsus$max_day))
print(paste("Humeral Length: Maximum growth rate is", max_growth_humeral_length$max_growth_rate, "on day", max_growth_humeral_length$max_day))
print(paste("Wing Length: Maximum growth rate is", max_growth_wing_length$max_growth_rate, "on day", max_growth_wing_length$max_day))

# Ensure the 'Day' column is numeric
growth_rate_summary_mass$Day <- as.numeric(growth_rate_summary_mass$Day)
growth_rate_summary_tarsus$Day <- as.numeric(growth_rate_summary_tarsus$Day)
growth_rate_summary_humeral_length$Day <- as.numeric(growth_rate_summary_humeral_length$Day)
growth_rate_summary_wing_length$Day <- as.numeric(growth_rate_summary_wing_length$Day)

# Add a variable identifier to each data frame
growth_rate_summary_mass$Variable <- "Mass"
growth_rate_summary_tarsus$Variable <- "Tarsus"
growth_rate_summary_humeral_length$Variable <- "Humeral Length"
growth_rate_summary_wing_length$Variable <- "Wing Length"

# Combine all data frames into one
combined_growth_rate_summary <- bind_rows(
  growth_rate_summary_mass,
  growth_rate_summary_tarsus,
  growth_rate_summary_humeral_length,
  growth_rate_summary_wing_length
)

# Identify the maximum growth rate and corresponding day for each variable
max_growth_rate_df <- combined_growth_rate_summary %>%
  group_by(Variable) %>%
  summarize(
    Max_Day = Day[which.max(Mean)],
    Max_Growth_Rate = max(Mean)
  ) %>%
  ungroup()

# Print the resulting data frame
print(max_growth_rate_df)

# Calculate mean growth rates for comaprison to other species in Table 2
# Function to calculate the mean growth rate across all time points
calculate_mean_growth_rate <- function(posterior_samples, time_points) {
  growth_rate_samples <- sapply(time_points, function(t) {
    apply(posterior_samples, 1, function(params) {
      gompertz_derivative(t, params["a"], params["b"], params["c"])
    })
  })
  
  # Calculate the mean growth rate over the entire study period (mean of all time points)
  mean_growth_rate <- rowMeans(growth_rate_samples)
  return(mean_growth_rate)
}

# Set the time points (for example, days 1 to 15)
time_points <- 1:15

# Calculate the mean growth rate for each measurement
mean_growth_rate_mass <- calculate_mean_growth_rate(posterior_samples_mass, time_points)
mean_growth_rate_tarsus <- calculate_mean_growth_rate(posterior_samples_tarsus, time_points)
mean_growth_rate_humeral_length <- calculate_mean_growth_rate(posterior_samples_humeral_length, time_points)
mean_growth_rate_wing_length <- calculate_mean_growth_rate(posterior_samples_wing_length, time_points)

# Summarize the results (mean and 95% credible intervals)
summarize_mean_growth_rate <- function(mean_growth_rate) {
  c(
    Mean = mean(mean_growth_rate),
    Lower_CI = quantile(mean_growth_rate, probs = 0.025),
    Upper_CI = quantile(mean_growth_rate, probs = 0.975)
  )
}

# Get the summarized mean growth rates for each measurement
mean_growth_rate_summary_mass <- summarize_mean_growth_rate(mean_growth_rate_mass)
mean_growth_rate_summary_tarsus <- summarize_mean_growth_rate(mean_growth_rate_tarsus)
mean_growth_rate_summary_humeral_length <- summarize_mean_growth_rate(mean_growth_rate_humeral_length)
mean_growth_rate_summary_wing_length <- summarize_mean_growth_rate(mean_growth_rate_wing_length)

# Print the results
print(mean_growth_rate_summary_mass)
print(mean_growth_rate_summary_tarsus)
print(mean_growth_rate_summary_humeral_length)
print(mean_growth_rate_summary_wing_length)

# Plot growth trajectories and growth rates for the top 4 variables

p <- ggplot(pred_data_all %>% filter(!is.na(Mean)),  
            aes(x = age_days, y = Mean, group = Measurement)) +
  geom_line(aes(linetype = Measurement), size = 1, color = "black") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "gray", alpha = 0.4) +
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

p_growth_rates <- ggplot(combined_growth_rate_summary %>% filter(!is.na(Mean)),  
                         aes(x = Day, y = Mean, group = Variable)) +
  geom_line(aes(linetype = Variable), size = 1, color = "black") +
  geom_ribbon(aes(ymin = Lower_95, ymax = Upper_95), fill = "gray", alpha = 0.4) +
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
  coord_cartesian(clip = "off")

# Combine the two plots into a 2-panel figure
combined_plot <- p + p_growth_rates + plot_layout(ncol = 1)

# Display the combined plot
print(combined_plot)

######################################################################################################################

# Generate the predictive table

# Assuming you have new_data with age_days and BandId (where BandId includes new levels)
new_data <- data.frame(age_days = 0:15, BandId = factor(rep("new_band", 16)))

# Generate predicted masses, allowing for new levels in BandId
predicted_variable <- posterior_epred(fit_mass, newdata = new_data, allow_new_levels = TRUE) # Replace fit_mass with other variables from the top 4

# Calculate the mean and 95% credible intervals for the predictions
variable_summary <- apply(predicted_variable, 2, function(x) {
  c(mean = mean(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
})

# Convert the summary statistics to a dataframe
variable_summary_df <- as.data.frame(t(variable_summary))

# Add age_days to the summary dataframe
variable_summary_df$age_days <- new_data$age_days

# Ensure the columns are in the desired order
variable_summary_df <- variable_summary_df %>%
  select(age_days, mean, `lower.2.5%`, `upper.97.5%`)

# Renaming the columns
names(variable_summary_df) <- c("Age_Days", "Predicted_Mean_Mass", "Lower_CI", "Upper_CI") # Replace Predicted_Mean_Mass with other variables from the top 4

# Viewing the modified table
print(variable_summary_df)

# Write the table to a CSV file
#write.csv(variable_summary_df, "ewpw_predictive_chick_mass.csv", row.names = FALSE) # Save a csv file for each variable, name however you like

# Combine tables

# Read in the csv files made above
pred_mass <- read.csv("ewpw_predictive_chick_mass.csv") # Use the csv files named above
pred_tarsus <- read.csv("ewpw_predictive_chick_tarsus.csv") # Use the csv files named above
pred_humeral_length <- read.csv("ewpw_predictive_chick_humeral_length.csv") # Use the csv files named above
pred_wing_length <- read.csv("ewpw_predictive_chick_wing_length.csv") # Use the csv files named above

# Join the data frames
combined_df <- pred_mass %>%
  inner_join(pred_tarsus, by = "Age_Days") %>%
  inner_join(pred_humeral_length, by = "Age_Days") %>%
  inner_join(pred_wing_length, by = "Age_Days")

combined_df <- combined_df %>%
  mutate(across(c(Predicted_Mean_Mass, Predicted_Mean_Tarsus, Predicted_Mean_Humeral_Length, Predicted_Mean_Wing_Length, 
                  Lower_CI.x, Upper_CI.x, Lower_CI.y, Upper_CI.y, Lower_CI.x.x, Upper_CI.x.x, Lower_CI.y.y, Upper_CI.y.y),~ round(., digits = 1)))

combined_df <- combined_df %>%
  mutate(Mass = paste0(formatC(Predicted_Mean_Mass, format = "f", digits = 1), 
                         " (", formatC(Lower_CI.x, format = "f", digits = 1), ", ", 
                         formatC(Upper_CI.x, format = "f", digits = 1), ")"))

combined_df <- combined_df %>%
  mutate(Tarsus = paste0(formatC(Predicted_Mean_Tarsus, format = "f", digits = 1), 
                         " (", formatC(Lower_CI.y, format = "f", digits = 1), ", ", 
                         formatC(Upper_CI.y, format = "f", digits = 1), ")"))

combined_df <- combined_df %>%
  mutate(HumeralLength = paste0(formatC(Predicted_Mean_Humeral_Length, format = "f", digits = 1), 
                                " (", formatC(Lower_CI.x.x, format = "f", digits = 1), ", ", 
                                formatC(Upper_CI.x.x, format = "f", digits = 1), ")"))

combined_df <- combined_df %>%
  mutate(WingLength = paste0(formatC(Predicted_Mean_Wing_Length, format = "f", digits = 1), 
                             " (", formatC(Lower_CI.y.y, format = "f", digits = 1), ", ", 
                             formatC(Upper_CI.y.y, format = "f", digits = 1), ")"))

combined_df <- combined_df %>%
  select(Age_Days, Mass, Tarsus, HumeralLength, WingLength)

# Write the data frame to a CSV file
#write.csv(combined_df, "comnined_chick_predictive_table.csv", row.names = FALSE) # Save and name however you like

# Predictive table combining all variables into a single prediction (adding all the posterior predicitons for the 4 variables together and then dividing by 4 to get a single unitless value)

new_data <- data.frame(age_days = 0:15, BandId = factor(rep("new_band", 16)))

# Generate predicted masses, allowing for new levels in BandId
predicted_variable_mass <- posterior_epred(fit_mass, newdata = new_data, allow_new_levels = TRUE)
predicted_variable_tarsus <- posterior_epred(fit_tarsus, newdata = new_data, allow_new_levels = TRUE)
predicted_variable_humeral_length <- posterior_epred(fit_humeral_length, newdata = new_data, allow_new_levels = TRUE)
predicted_variable_wing_length <- posterior_epred(fit_wing_length, newdata = new_data, allow_new_levels = TRUE)

# Sum the predictions across the four models and average them
combined_predictions <- (predicted_variable_mass + 
                           predicted_variable_tarsus + 
                           predicted_variable_humeral_length + 
                           predicted_variable_wing_length) / 4

# Step 1: Calculate the mean and 95% credible intervals for each day (column) across all posterior draws
summary_stats <- apply(combined_predictions, 2, function(column) {
  c(
    Mean = mean(column),
    Lower_CI = quantile(column, probs = 0.025),
    Upper_CI = quantile(column, probs = 0.975)
  )
})

# Convert the summary statistics to a data frame
summary_stats_df <- as.data.frame(summary_stats)

# Transpose the data frame so that each row corresponds to a day
summary_stats_df <- t(summary_stats_df)

# Add age_days to the summary dataframe
summary_stats_df <- data.frame(Age_Days = 0:15, summary_stats_df)

# View the resulting data frame
print(summary_stats_df)

# Write the data frame to a CSV file
#write.csv(summary_stats_df, "comnined_chick_predictive_table_all_variables_combined.csv", row.names = FALSE) # Save and name however you like

######################################################################################################################

# Generate posterior predictions for the other variables (CV greater than 15%): bill_length_from_nares, posterior_bill_width, anterior_bill_width, mouth_width, tail_length

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

# Generate a new data frame for the prediction range
new_data <- data.frame(age_days = age_seq)

# Use the 'posterior_predict' function to get draws from the posterior predictive distribution
posterior_draws <- posterior_predict(fit_tail_length, newdata = new_data, re_formula = NA)

# Calculate the mean and credible interval for the predictions
pred_summary <- apply(posterior_draws, 2, function(x) {
  c(mean = mean(x), ci_lower = quantile(x, probs = 0.025), ci_upper = quantile(x, probs = 0.975))
})

# Convert summary to data frame for plotting
pred_summary <- t(pred_summary) 
pred_data <- as.data.frame(pred_summary)  
names(pred_data) <- c("mean", "ci_lower", "ci_upper")  
pred_data$age_days <- age_seq 

# Create prediction data for each measurement

# Bill Length from Nares
posterior_draws_bill_length_from_nares <- posterior_predict(fit_bill_length_from_nares, newdata = new_data, re_formula = NA)
pred_summary_bill_length_from_nares <- apply(posterior_draws_bill_length_from_nares, 2, function(x) {
  c(mean = mean(x), ci_lower = quantile(x, probs = 0.025), ci_upper = quantile(x, probs = 0.975))
})
pred_data_bill_length_from_nares <- as.data.frame(t(pred_summary_bill_length_from_nares))
pred_data_bill_length_from_nares$age_days <- age_seq
pred_data_bill_length_from_nares$Measurement <- "Bill Length from Nares"

# Posterior Bill Width
posterior_draws_posterior_bill_width <- posterior_predict(fit_posterior_bill_width, newdata = new_data, re_formula = NA)
pred_summary_posterior_bill_width <- apply(posterior_draws_posterior_bill_width, 2, function(x) {
  c(mean = mean(x), ci_lower = quantile(x, probs = 0.025), ci_upper = quantile(x, probs = 0.975))
})
pred_data_posterior_bill_width <- as.data.frame(t(pred_summary_posterior_bill_width))
pred_data_posterior_bill_width$age_days <- age_seq
pred_data_posterior_bill_width$Measurement <- "Posterior Bill Width"

# Anterior Bill Width
posterior_draws_anterior_bill_width <- posterior_predict(fit_anterior_bill_width, newdata = new_data, re_formula = NA)
pred_summary_anterior_bill_width <- apply(posterior_draws_anterior_bill_width, 2, function(x) {
  c(mean = mean(x), ci_lower = quantile(x, probs = 0.025), ci_upper = quantile(x, probs = 0.975))
})
pred_data_anterior_bill_width <- as.data.frame(t(pred_summary_anterior_bill_width))
pred_data_anterior_bill_width$age_days <- age_seq
pred_data_anterior_bill_width$Measurement <- "Anterior Bill Width"

# Mouth Width
posterior_draws_mouth_width <- posterior_predict(fit_mouth_width, newdata = new_data, re_formula = NA)
pred_summary_mouth_width <- apply(posterior_draws_mouth_width, 2, function(x) {
  c(mean = mean(x), ci_lower = quantile(x, probs = 0.025), ci_upper = quantile(x, probs = 0.975))
})
pred_data_mouth_width <- as.data.frame(t(pred_summary_mouth_width))
pred_data_mouth_width$age_days <- age_seq
pred_data_mouth_width$Measurement <- "Mouth Width"

# Tail Length
posterior_draws_tail_length <- posterior_predict(fit_tail_length, newdata = new_data, re_formula = NA)
pred_summary_tail_length <- apply(posterior_draws_tail_length, 2, function(x) {
  c(mean = mean(x), ci_lower = quantile(x, probs = 0.025), ci_upper = quantile(x, probs = 0.975))
})
pred_data_tail_length <- as.data.frame(t(pred_summary_tail_length))
pred_data_tail_length$age_days <- age_seq
pred_data_tail_length$Measurement <- "Tail Length"

# Combine all predictions into one dataframe
combined_pred_data <- rbind(pred_data_bill_length_from_nares, pred_data_posterior_bill_width, pred_data_anterior_bill_width, pred_data_mouth_width, pred_data_tail_length)

# Rename columns to ensure consistency
colnames(combined_pred_data) <- c("mean", "ci_lower", "ci_upper", "age_days", "Measurement")

# Reshape other_data to long format
observed_data_long <- other_data %>%
  pivot_longer(
    cols = c(culmen, bill_width_bn, bill_width_an, mouth_width, tail_length),  # Columns to reshape
    names_to = "Measurement",  # New column to hold measurement names
    values_to = "value"  # New column to hold measurement values
  )

combined_5_measurements <- ggplot(combined_pred_data, aes(x = age_days, y = mean, shape = Measurement)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "gray80") +  # Prediction intervals
  geom_line(size = 1, colour = "gray20", linetype = "solid") +  # Fitted growth curves with solid lines
  geom_point(data = observed_data_long, aes(x = age_days, y = value, shape = Measurement), colour = "black", alpha = 0.7, size = 2.5) +  # Observed data points with different shapes
  labs(
    x = "Days",
    y = "Millimeters",
    shape = "Measurement"  # Add a label for the legend
  ) +
  facet_wrap(~ Measurement, scales = "free_y") +  # Facet by measurement type
  scale_shape_manual(values = c("Bill Length from Nares" = 16, "Posterior Bill Width" = 17, "Anterior Bill Width" = 15, "Mouth Width" = 18, "Tail Length" = 3)) +  # Custom shapes for each measurement
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 26, margin = margin(t = 12)),
    axis.title.y = element_text(size = 26, margin = margin(r = 12)),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    strip.text = element_blank(),  # Remove facet plot titles
    axis.line = element_line(size = 0.5),
    legend.position = "right",  # Position the legend on the right
    legend.text = element_text(size = 18),  # Increase legend text size
    legend.title = element_text(size = 20)  # Increase legend title size
  )

print(combined_5_measurements)
