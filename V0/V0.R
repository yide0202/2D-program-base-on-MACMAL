

###### Version 0.6 Cluster Survival 


# Assume we have a vector `death_counts` representing the number of deaths for each time interval
# `total_pop` represents the population at risk or the total population for each interval
# We aim to fit a segmented model to capture the heterogeneity in mortality rates

# First, install Rcpp if not already installed
# install.packages("Rcpp")

library(Rcpp)

# Embed C++ code in R
sourceCpp("cluster_survival.cpp")

# Example data (for demonstration purposes)
set.seed(123)
N <- 100
# Generate random death counts with varying mortality rates
death_counts <- rpois(N, lambda = c(rep(1,40), rep(5,20), rep(1,40)))
total_pop <- rep(1000, N)  # Population at risk is constant for simplicity

# Run the C++ function for model fitting and clustering analysis
result <- runClusterSurvival(death_counts, total_pop)

# `result` contains the selected model, estimated rates for each interval,
# interval boundaries, and model-averaged results
print(result)

# Compute the observed mortality rate for each interval
obs_rate <- death_counts / total_pop

# Plot the observed mortality rates as a scatter plot
plot(obs_rate, type = "p", pch = 16, col = "gray50",
     xlab = "Interval Index", ylab = "Mortality Rate", 
     main = "Mortality Rate and Identified Heterogeneous Intervals")

# Add a smoothing line to observe the overall trend
lines(smooth.spline(obs_rate), col = "blue", lwd = 2)

# Highlight the detected intervals using the results
cs <- result$cs  # Start indices of clusters
ce <- result$ce  # End indices of clusters
p0 <- result$p0  # External rates
pc <- result$pc  # Cluster-specific rates

# Mark the identified clusters on the plot
for (i in seq_along(cs)) {
  start_pos <- cs[i]
  end_pos <- ce[i]
  
  # Draw shaded rectangles to highlight clusters
  rect(xleft = start_pos, ybottom = min(obs_rate) * 0.9, 
       xright = end_pos, ytop = max(obs_rate) * 1.1, 
       border = NA, col = rgb(1, 0, 0, alpha = 0.1))
  
  # Add horizontal lines representing `pc` within clusters
  segments(x0 = start_pos, y0 = pc[i], x1 = end_pos, y1 = pc[i],
           col = "red", lwd = 2)
  
  # Optionally add reference lines for `p0` (external rates)
  segments(x0 = start_pos, y0 = p0[i], x1 = end_pos, y1 = p0[i],
           col = "green", lwd = 2, lty = 2)
}

# Add a legend to explain the plot
legend("topright", legend = c("Observed Mortality Rate", "Smoothed Trend Line", 
                              "Detected Clusters", "Cluster Rate (pc)", "External Rate (p0)"),
       col = c("gray50", "blue", "red", "red", "green"), 
       pch = c(16, NA, NA, NA, NA), 
       lty = c(NA, 1, NA, 1, 2), 
       lwd = c(1, 2, NA, 2, 2), 
       bty = "n")

# Assign risk scores based on the ratio of `pc` to `p0`
ratio <- pc / p0
scaling_factor <- 2  # Adjust based on data distribution
scores <- 5 + scaling_factor * log2(ratio)
scores[scores > 10] <- 10  # Cap the score at 10
scores[scores < 1] <- 1  # Set a minimum score of 1

# Annotate the plot with scores
for (i in seq_along(cs)) {
  start_pos <- cs[i]
  end_pos <- ce[i]
  mid_pos <- (start_pos + end_pos) / 2
  
  # Add score annotations to the plot
  text(x = mid_pos, y = pc[i] + (max(obs_rate) * 0.05),
       labels = paste0("Score: ", round(scores[i], 1)), 
       col = "red", font = 2)
}

# Inspect the structure of the result object
str(result)




# Load required libraries
library(ggplot2)
library(dplyr)

# Generate example data
set.seed(123)
N <- 100
death_counts <- rpois(N, lambda = c(rep(1, 40), rep(5, 20), rep(1, 40)))  # Simulated death counts
total_pop <- rep(1000, N)  # Population at risk
result <- runClusterSurvival(death_counts, total_pop)  # Run clustering function

# Extract results from the clustering output
cs <- result$cs  # Start positions of clusters
ce <- result$ce  # End positions of clusters
p0 <- result$p0  # Background death rate outside clusters
pc <- result$pc  # Death rate within clusters

# Create a data frame for visualization
data <- data.frame(
  Interval = 1:N,
  DeathRate = death_counts / total_pop,  # Observed death rate
  Deaths = death_counts  # Death counts
)

# Mark clustered intervals
data$Cluster <- 0  # Default: no cluster
for (i in seq_along(cs)) {
  data$Cluster[data$Interval >= cs[i] & data$Interval <= ce[i]] <- i
}

# Separate clustered and non-clustered data
clusters <- data %>% filter(Cluster > 0)
non_clusters <- data %>% filter(Cluster == 0)

# ============================
# Plot 1: Death Rate with Cluster Annotations
# ============================
p1 <- ggplot() +
  geom_point(data = non_clusters, aes(x = Interval, y = DeathRate), 
             color = "gray", alpha = 0.7, size = 2) +  # Non-clustered points
  geom_point(data = clusters, aes(x = Interval, y = DeathRate, color = factor(Cluster)), 
             size = 3, alpha = 0.8) +  # Highlighted clustered points
  geom_rect(data = data.frame(cs, ce), 
            aes(xmin = cs, xmax = ce, ymin = 0, ymax = max(data$DeathRate) * 1.1),
            fill = "red", alpha = 0.1, inherit.aes = FALSE) +  # Highlight cluster regions
  geom_text(data = data.frame(cs, ce), 
            aes(x = cs, y = max(data$DeathRate) * 1.05, label = paste("Start:", cs)), 
            color = "red", angle = 90, hjust = 1) +  # Start annotation
  geom_text(data = data.frame(cs, ce), 
            aes(x = ce, y = max(data$DeathRate) * 1.05, label = paste("End:", ce)), 
            color = "red", angle = 90, hjust = 0) +  # End annotation
  annotate("text", x = 25, y = 0.006, 
           label = "Highlighted clusters indicate regions\nwith significantly different death rates", 
           color = "darkred", size = 4, hjust = 0) +
  labs(title = "Observed Death Rate with Cluster Annotations",
       subtitle = "Highlighted regions represent detected clusters",
       x = "Time Interval", y = "Death Rate") +
  scale_color_brewer(palette = "Set1", name = "Cluster ID") +
  theme_minimal()

# ============================
# Plot 2: Scatterplot of Death Counts by Cluster
# ============================
p2 <- ggplot(data, aes(x = Interval, y = Deaths)) +
  geom_point(aes(size = Deaths, color = factor(Cluster)), alpha = 0.8) +  # Scatterplot with death counts
  geom_text(aes(label = ifelse(Deaths > 5, Deaths, "")), 
            color = "black", size = 3, vjust = -1) +  # Annotate high death counts
  scale_size_continuous(range = c(1, 5), name = "Deaths") +
  scale_color_manual(values = c("gray", rainbow(length(unique(data$Cluster)) - 1)),
                     name = "Cluster ID") +
  labs(title = "Scatterplot of Death Events by Cluster",
       subtitle = "Larger points indicate higher death counts",
       x = "Time Interval", y = "Death Counts") +
  theme_minimal()

# ============================
# Plot 3: Smoothed Death Rate Trend with Cluster Highlights
# ============================


p3 <- ggplot(data, aes(x = Interval, y = DeathRate)) +
  geom_point(aes(color = DeathRate), size = 2, alpha = 0.8) +  # Observed points
  geom_line(aes(y = RollingAvg), color = "blue", size = 1.2) +  # Rolling average line
  geom_rect(data = data.frame(cs, ce), 
            aes(xmin = cs, xmax = ce, ymin = 0, ymax = max(data$DeathRate) * 1.1),
            fill = "red", alpha = 0.1, inherit.aes = FALSE) +  # Highlight clusters
  scale_color_gradient(low = "blue", high = "red", name = "Death Rate") +
  annotate("text", x = 50, y = max(data$DeathRate) * 1.05, 
           label = "Rolling average of death rate trend\nwith highlighted clusters", 
           color = "blue", size = 4, hjust = 0.5) +
  labs(title = "Death Rate Trend Using Rolling Average",
       x = "Time Interval", y = "Death Rate") +
  theme_minimal()


library(zoo)

# Calculate rolling average with a window size of 5
data$RollingAvg <- rollmean(data$DeathRate, k = 5, fill = NA)

# Update plot with rolling average
p3 <- ggplot(data, aes(x = Interval, y = DeathRate)) +
  geom_point(aes(color = DeathRate), size = 2, alpha = 0.8) +  # Observed points
  geom_line(aes(y = RollingAvg), color = "blue", size = 1.2) +  # Rolling average line
  geom_rect(data = data.frame(cs, ce), 
            aes(xmin = cs, xmax = ce, ymin = 0, ymax = max(data$DeathRate) * 1.1),
            fill = "red", alpha = 0.1, inherit.aes = FALSE) +  # Highlight clusters
  scale_color_gradient(low = "blue", high = "red", name = "Death Rate") +
  annotate("text", x = 50, y = max(data$DeathRate) * 1.05, 
           label = "Rolling average of death rate trend\nwith highlighted clusters", 
           color = "blue", size = 4, hjust = 0.5) +
  labs(title = "Death Rate Trend Using Rolling Average",
       x = "Time Interval", y = "Death Rate") +
  theme_minimal()







# ============================
# Output: Display each plot separately
# ============================
print(p1)  # Plot 1: Death Rate with Cluster Annotations
print(p2)  # Plot 2: Scatterplot of Death Counts
print(p3)  # Plot 3: Smoothed Death Rate Trend with Clusters



p3 <- ggplot(data, aes(x = Interval, y = DeathRate)) +
  geom_point(aes(color = DeathRate), size = 2, alpha = 0.8) +  # Observed points
  geom_smooth(se = FALSE, method = "loess", span = 0.2, color = "blue", size = 1.2) + # Adjust span for peaks
  geom_rect(data = data.frame(cs, ce), 
            aes(xmin = cs, xmax = ce, ymin = 0, ymax = max(data$DeathRate) * 1.1),
            fill = "red", alpha = 0.1, inherit.aes = FALSE) +  # Highlight clusters
  scale_color_gradient(low = "blue", high = "red", name = "Death Rate") +
  annotate("text", x = 50, y = max(data$DeathRate) * 1.05, 
           label = "Smoothed death rate trend\nwith highlighted clusters", 
           color = "blue", size = 4, hjust = 0.5) +
  labs(title = "Improved Death Rate Trend with Localized Peaks",
       x = "Time Interval", y = "Death Rate") +
  theme_minimal()



library(zoo)

# Calculate rolling average with a window size of 5
data$RollingAvg <- rollmean(data$DeathRate, k = 5, fill = NA)

# Update plot with rolling average
p3 <- ggplot(data, aes(x = Interval, y = DeathRate)) +
  geom_point(aes(color = DeathRate), size = 2, alpha = 0.8) +  # Observed points
  geom_line(aes(y = RollingAvg), color = "blue", size = 1.2) +  # Rolling average line
  geom_rect(data = data.frame(cs, ce), 
            aes(xmin = cs, xmax = ce, ymin = 0, ymax = max(data$DeathRate) * 1.1),
            fill = "red", alpha = 0.1, inherit.aes = FALSE) +  # Highlight clusters
  scale_color_gradient(low = "blue", high = "red", name = "Death Rate") +
  annotate("text", x = 50, y = max(data$DeathRate) * 1.05, 
           label = "Rolling average of death rate trend\nwith highlighted clusters", 
           color = "blue", size = 4, hjust = 0.5) +
  labs(title = "Death Rate Trend Using Rolling Average",
       x = "Time Interval", y = "Death Rate") +
  theme_minimal()








