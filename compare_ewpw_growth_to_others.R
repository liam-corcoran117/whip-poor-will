# Load necessary libraries
library(nlstools)
library(ggplot2)

# Data for Fiery-necked Nightjar
nightjar_data <- data.frame(
  day = 1:20,
  mass = c(5.3, 6.2, 8.2, 10.2, 13.1, 15.4, 18.6, 20.9, 23.4, 25.5, 27.2, 29.4, 30.1, 33.0, 34.7, 35.2, 37.1, 36.6, 39.1, 40.0),
  tarsus = c(7.5, 10.2, 10.1, 12.1, 12.2, 12.6, 13.4, 13.6, 13.3, 15.1, 15.0, 16.5, 15.6, 15.3, 15.9, 17.5, 17.0, 17.3, 19.0, 18.0),
  wing = c(10.3, 13.6, 14.0, 16.9, 21.9, 28.3, 35.3, 41.1, 47.9, 53.1, 61.3, 66.5, 71.7, 77.7, 83.6, 88.1, 93.7, 94.4, 95.3, 102.0)
)

# Data for Eastern Kingbird
eastern_kingbird_data <- data.frame(
  day = 1:14,
  mass = c(3.58, 5.53, 7.82, 10.36, 13.54, 16.56, 19.39, 22.27, 25.21, 26.61, 27.59, 28.80, 29.77, 29.79),
  tarsus = c(6.43, 7.52, 8.74, 10.24, 11.86, 13.29, 14.61, 15.72, 16.59, 17.16, 17.61, 17.91, 18.17, 18.26),
  wing = c(0, 0, 0, 0.86, 2.6, 5.19, 8.49, 12.43, 16.34, 20.43, 24.36, 28.15, 31.89, 35.13)  # Wing length
)

# Data for Eastern Phoebe
eastern_phoebe_data <- data.frame(
  day = 1:14,
  mass = c(1.91, 2.77, 3.98, 5.43, 7.10, 8.89, 10.96, 12.51, 14.22, 15.36, 16.32, 16.91, 17.16, 17.47),
  tarsus = c(5.47, 6.25, 7.18, 8.46, 9.78, 11.26, 12.63, 13.82, 14.91, 15.86, 16.43, 16.93, 17.26, 17.41),
  wing = c(0, 0, 0, 0.43, 1.95, 4.24, 7.30, 10.70, 14.27, 18.13, 21.87, 25.41, 28.60, 31.59)  # Wing length
)

# Define Gompertz model
gompertz <- function(t, a, b, c) {
  a * exp(-b * exp(-c * t))
}

# Fit Gompertz models for each species and each measurement
# Priors for parameter a derived from adult size values for each species found in the Birds of the World accounts
# Fiery-necked Nightjar
fit_mass_nightjar <- nls(mass ~ gompertz(day, a, b, c), data = nightjar_data, start = list(a = 50, b = 1, c = 0.1))
fit_tarsus_nightjar <- nls(tarsus ~ gompertz(day, a, b, c), data = nightjar_data, start = list(a = 14, b = 1, c = 0.1))
fit_wing_nightjar <- nls(wing ~ gompertz(day, a, b, c), data = nightjar_data, start = list(a = 155, b = 1, c = 0.1))

# Eastern Kingbird
fit_mass_kingbird <- nls(mass ~ gompertz(day, a, b, c), data = eastern_kingbird_data, start = list(a = 40, b = 1, c = 0.1))
fit_tarsus_kingbird <- nls(tarsus ~ gompertz(day, a, b, c), data = eastern_kingbird_data, start = list(a = 19, b = 1, c = 0.1))
fit_wing_kingbird <- nls(wing ~ gompertz(day, a, b, c), data = eastern_kingbird_data, start = list(a = 115, b = 1, c = 0.1))

# Eastern Phoebe
fit_mass_phoebe <- nls(mass ~ gompertz(day, a, b, c), data = eastern_phoebe_data, start = list(a = 20, b = 1, c = 0.1))
fit_tarsus_phoebe <- nls(tarsus ~ gompertz(day, a, b, c), data = eastern_phoebe_data, start = list(a = 18, b = 1, c = 0.1))
fit_wing_phoebe <- nls(wing ~ gompertz(day, a, b, c), data = eastern_phoebe_data, start = list(a = 100, b = 1, c = 0.1))

# Extract parameters for each model
params_mass_nightjar <- coef(fit_mass_nightjar)
params_tarsus_nightjar <- coef(fit_tarsus_nightjar)
params_wing_nightjar <- coef(fit_wing_nightjar)

params_mass_kingbird <- coef(fit_mass_kingbird)
params_tarsus_kingbird <- coef(fit_tarsus_kingbird)
params_wing_kingbird <- coef(fit_wing_kingbird)

params_mass_phoebe <- coef(fit_mass_phoebe)
params_tarsus_phoebe <- coef(fit_tarsus_phoebe)
params_wing_phoebe <- coef(fit_wing_phoebe)

# Define the first derivative of the Gompertz equation (for growth rate)
gompertz_derivative <- function(t, a, b, c) {
  a * b * c * exp(-b * exp(-c * t)) * exp(-c * t)
}

# Function to calculate growth rate summary
calculate_growth_rate_summary <- function(params, time_points) {
  growth_rate_samples <- sapply(time_points, function(t) {
    gompertz_derivative(t, params["a"], params["b"], params["c"])
  })
  
  data.frame(Day = time_points, GrowthRate = growth_rate_samples)
}

# Calculate growth rate summaries for all measurements
time_points_nightjar <- 1:20
growth_rate_summary_mass_nightjar <- calculate_growth_rate_summary(params_mass_nightjar, time_points_nightjar)
growth_rate_summary_tarsus_nightjar <- calculate_growth_rate_summary(params_tarsus_nightjar, time_points_nightjar)
growth_rate_summary_wing_nightjar <- calculate_growth_rate_summary(params_wing_nightjar, time_points_nightjar)

time_points_kingbird_phoebe <- 1:14
growth_rate_summary_mass_kingbird <- calculate_growth_rate_summary(params_mass_kingbird, time_points_kingbird_phoebe)
growth_rate_summary_tarsus_kingbird <- calculate_growth_rate_summary(params_tarsus_kingbird, time_points_kingbird_phoebe)
growth_rate_summary_wing_kingbird <- calculate_growth_rate_summary(params_wing_kingbird, time_points_kingbird_phoebe)

growth_rate_summary_mass_phoebe <- calculate_growth_rate_summary(params_mass_phoebe, time_points_kingbird_phoebe)
growth_rate_summary_tarsus_phoebe <- calculate_growth_rate_summary(params_tarsus_phoebe, time_points_kingbird_phoebe)
growth_rate_summary_wing_phoebe <- calculate_growth_rate_summary(params_wing_phoebe, time_points_kingbird_phoebe)

# Combine all growth rate summaries for plotting
growth_rate_summary_mass_kingbird$Variable <- "Kingbird Mass"
growth_rate_summary_tarsus_kingbird$Variable <- "Kingbird Tarsus"
growth_rate_summary_wing_kingbird$Variable <- "Kingbird Wing"

growth_rate_summary_mass_phoebe$Variable <- "Phoebe Mass"
growth_rate_summary_tarsus_phoebe$Variable <- "Phoebe Tarsus"
growth_rate_summary_wing_phoebe$Variable <- "Phoebe Wing"

growth_rate_summary_mass_nightjar$Variable <- "Nightjar Mass"
growth_rate_summary_tarsus_nightjar$Variable <- "Nightjar Tarsus"
growth_rate_summary_wing_nightjar$Variable <- "Nightjar Wing"

combined_growth_rate_summary <- rbind(
  growth_rate_summary_mass_kingbird,
  growth_rate_summary_tarsus_kingbird,
  growth_rate_summary_wing_kingbird,
  growth_rate_summary_mass_phoebe,
  growth_rate_summary_tarsus_phoebe,
  growth_rate_summary_wing_phoebe,
  growth_rate_summary_mass_nightjar,
  growth_rate_summary_tarsus_nightjar,
  growth_rate_summary_wing_nightjar
)

# Function to find the maximum growth rates and corresponding day
find_max_growth_rate <- function(growth_rate_summary) {
  max_index <- which.max(growth_rate_summary$GrowthRate)
  max_day <- growth_rate_summary$Day[max_index]
  max_growth_rate <- growth_rate_summary$GrowthRate[max_index]
  list(max_day = max_day, max_growth_rate = max_growth_rate)
}

# Find the maximum growth rates and corresponding days for each variable
max_growth_mass_kingbird <- find_max_growth_rate(growth_rate_summary_mass_kingbird)
max_growth_tarsus_kingbird <- find_max_growth_rate(growth_rate_summary_tarsus_kingbird)
max_growth_wing_kingbird <- find_max_growth_rate(growth_rate_summary_wing_kingbird)

max_growth_mass_phoebe <- find_max_growth_rate(growth_rate_summary_mass_phoebe)
max_growth_tarsus_phoebe <- find_max_growth_rate(growth_rate_summary_tarsus_phoebe)
max_growth_wing_phoebe <- find_max_growth_rate(growth_rate_summary_wing_phoebe)

max_growth_mass_nightjar <- find_max_growth_rate(growth_rate_summary_mass_nightjar)
max_growth_tarsus_nightjar <- find_max_growth_rate(growth_rate_summary_tarsus_nightjar)
max_growth_wing_nightjar <- find_max_growth_rate(growth_rate_summary_wing_nightjar)

# Print the results
print(paste("Kingbird Mass: Maximum growth rate is", round(max_growth_mass_kingbird$max_growth_rate, 2), "on day", max_growth_mass_kingbird$max_day))
print(paste("Kingbird Tarsus: Maximum growth rate is", round(max_growth_tarsus_kingbird$max_growth_rate, 2), "on day", max_growth_tarsus_kingbird$max_day))
print(paste("Kingbird Wing: Maximum growth rate is", round(max_growth_wing_kingbird$max_growth_rate, 2), "on day", max_growth_wing_kingbird$max_day))

print(paste("Phoebe Mass: Maximum growth rate is", round(max_growth_mass_phoebe$max_growth_rate, 2), "on day", max_growth_mass_phoebe$max_day))
print(paste("Phoebe Tarsus: Maximum growth rate is", round(max_growth_tarsus_phoebe$max_growth_rate, 2), "on day", max_growth_tarsus_phoebe$max_day))
print(paste("Phoebe Wing: Maximum growth rate is", round(max_growth_wing_phoebe$max_growth_rate, 2), "on day", max_growth_wing_phoebe$max_day))

print(paste("Nightjar Mass: Maximum growth rate is", round(max_growth_mass_nightjar$max_growth_rate, 2), "on day", max_growth_mass_nightjar$max_day))
print(paste("Nightjar Tarsus: Maximum growth rate is", round(max_growth_tarsus_nightjar$max_growth_rate, 2), "on day", max_growth_tarsus_nightjar$max_day))
print(paste("Nightjar Wing: Maximum growth rate is", round(max_growth_wing_nightjar$max_growth_rate, 2), "on day", max_growth_wing_nightjar$max_day))

# Function to calculate the mean growth rate over the study period
calculate_mean_growth_rate <- function(growth_rate_summary) {
  mean(growth_rate_summary$GrowthRate)
}

# Calculate the mean growth rates for each variable
mean_growth_mass_kingbird <- calculate_mean_growth_rate(growth_rate_summary_mass_kingbird)
mean_growth_tarsus_kingbird <- calculate_mean_growth_rate(growth_rate_summary_tarsus_kingbird)
mean_growth_wing_kingbird <- calculate_mean_growth_rate(growth_rate_summary_wing_kingbird)

mean_growth_mass_phoebe <- calculate_mean_growth_rate(growth_rate_summary_mass_phoebe)
mean_growth_tarsus_phoebe <- calculate_mean_growth_rate(growth_rate_summary_tarsus_phoebe)
mean_growth_wing_phoebe <- calculate_mean_growth_rate(growth_rate_summary_wing_phoebe)

mean_growth_mass_nightjar <- calculate_mean_growth_rate(growth_rate_summary_mass_nightjar)
mean_growth_tarsus_nightjar <- calculate_mean_growth_rate(growth_rate_summary_tarsus_nightjar)
mean_growth_wing_nightjar <- calculate_mean_growth_rate(growth_rate_summary_wing_nightjar)

# Print the mean growth rates
print(paste("Kingbird Mass: Mean growth rate is", round(mean_growth_mass_kingbird, 2)))
print(paste("Kingbird Tarsus: Mean growth rate is", round(mean_growth_tarsus_kingbird, 2)))
print(paste("Kingbird Wing: Mean growth rate is", round(mean_growth_wing_kingbird, 2)))

print(paste("Phoebe Mass: Mean growth rate is", round(mean_growth_mass_phoebe, 2)))
print(paste("Phoebe Tarsus: Mean growth rate is", round(mean_growth_tarsus_phoebe, 2)))
print(paste("Phoebe Wing: Mean growth rate is", round(mean_growth_wing_phoebe, 2)))

print(paste("Nightjar Mass: Mean growth rate is", round(mean_growth_mass_nightjar, 2)))
print(paste("Nightjar Tarsus: Mean growth rate is", round(mean_growth_tarsus_nightjar, 2)))
print(paste("Nightjar Wing: Mean growth rate is", round(mean_growth_wing_nightjar, 2)))


# Plot for Fiery-necked Nightjar
nightjar_plot <- ggplot() +
  geom_line(data = growth_rate_summary_mass_nightjar, aes(x = Day, y = GrowthRate, color = "Nightjar Mass"), size = 1) +
  geom_line(data = growth_rate_summary_tarsus_nightjar, aes(x = Day, y = GrowthRate, color = "Nightjar Tarsus"), size = 1) +
  geom_line(data = growth_rate_summary_wing_nightjar, aes(x = Day, y = GrowthRate, color = "Nightjar Wing"), size = 1) +
  labs(title = "Growth Rate Over Time for Fiery-necked Nightjar",
       y = "Growth Rate (g/day or mm/day)", x = "Day") +
  theme_minimal() +
  scale_color_manual(values = c("Nightjar Mass" = "black", "Nightjar Tarsus" = "blue", "Nightjar Wing" = "red")) +
  theme(legend.title = element_blank())

print(nightjar_plot)

# Plot for Eastern Kingbird
kingbird_plot <- ggplot() +
  geom_line(data = growth_rate_summary_mass_kingbird, aes(x = Day, y = GrowthRate, color = "Kingbird Mass"), size = 1) +
  geom_line(data = growth_rate_summary_tarsus_kingbird, aes(x = Day, y = GrowthRate, color = "Kingbird Tarsus"), size = 1) +
  geom_line(data = growth_rate_summary_wing_kingbird, aes(x = Day, y = GrowthRate, color = "Kingbird Wing"), size = 1) +
  labs(title = "Growth Rate Over Time for Eastern Kingbird",
       y = "Growth Rate (g/day or mm/day)", x = "Day") +
  theme_minimal() +
  scale_color_manual(values = c("Kingbird Mass" = "black", "Kingbird Tarsus" = "blue", "Kingbird Wing" = "red")) +
  theme(legend.title = element_blank())

print(kingbird_plot)

# Plot for Eastern Phoebe
phoebe_plot <- ggplot() +
  geom_line(data = growth_rate_summary_mass_phoebe, aes(x = Day, y = GrowthRate, color = "Phoebe Mass"), size = 1) +
  geom_line(data = growth_rate_summary_tarsus_phoebe, aes(x = Day, y = GrowthRate, color = "Phoebe Tarsus"), size = 1) +
  geom_line(data = growth_rate_summary_wing_phoebe, aes(x = Day, y = GrowthRate, color = "Phoebe Wing"), size = 1) +
  labs(title = "Growth Rate Over Time for Eastern Phoebe",
       y = "Growth Rate (g/day or mm/day)", x = "Day") +
  theme_minimal() +
  scale_color_manual(values = c("Phoebe Mass" = "black", "Phoebe Tarsus" = "blue", "Phoebe Wing" = "red")) +
  theme(legend.title = element_blank())

print(phoebe_plot)
