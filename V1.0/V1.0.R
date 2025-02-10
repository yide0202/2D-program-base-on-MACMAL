library(Rcpp)
sourceCpp("cluster_survival_covariates.cpp")
# Generate small-scale data
N <- 10
set.seed(123)
pop_r    <- rep(50, N)
base_p   <- rep(0.05, N); base_p[4:6] <- 0.15  # Slightly higher in the middle period
deaths_r <- rbinom(N, pop_r, base_p)
age      <- round(seq(30,60,length.out=N))
gender   <- rbinom(N,1,0.4)
X        <- cbind(age, gender)

# Call the function
res <- runClusterSurvivalCovariates(deaths_r, pop_r, X, 0) # 0=BIC
print(res)


library(ggplot2)

# 1. Generate simulated data
N <- 15
set.seed(202309)   # Set random seed for reproducibility

# At-risk population
pop_r <- rep(50, N)

# Two covariates: age gradually increasing, gender random 0/1
age    <- round(seq(30, 60, length.out=N))
gender <- rbinom(N, 1, 0.4)

# Construct covariate matrix X
X <- cbind(age, gender)

# Manually set "true" mortality rate: higher in the interval (6~10), lower elsewhere
true_p <- rep(0.05, N)   # baseline
true_p[6:10] <- 0.15     # Slightly higher in the middle range

# Generate number of deaths
deaths_r <- rbinom(N, size=pop_r, prob=true_p)

# Create a data frame
df <- data.frame(
  Time     = 1:N,
  Deaths   = deaths_r,
  Pop      = pop_r,
  Observed = deaths_r / pop_r,
  Age      = age,
  Gender   = gender,
  TrueRate = true_p
)
print(df)

# 2. Call runClusterSurvivalCovariates (from cluster_survival_covariates.cpp)
#    Ensure the function is successfully exported to the R environment
res <- runClusterSurvivalCovariates(deaths_r, pop_r, X, criterion_type=0) 
# criterion_type=0 => BIC, 1 => AIC

# View results
print(res)

# Extract cluster range
cs <- res$cs
ce <- res$ce
cat("Detected cluster range: [", cs, ", ", ce, "]\n")

# 3. Visualization
# Using ggplot2 to show:
#   - x-axis: Time
#   - y-axis: Observed mortality rate
#   - Light red rectangle: cluster range detected by the model (cs ~ ce)
#   - Overlay true mortality rate TrueRate with lines and points

ggplot(df, aes(x=Time)) +
  # Observed mortality rate -> bar or dot plot
  geom_col(aes(y=Observed), fill="steelblue", alpha=0.6) +
  
  # Highlight cluster range on the plot (light red background)
  geom_rect(
    aes(xmin=cs-0.4, xmax=ce+0.4, ymin=-Inf, ymax=Inf),
    fill="red", alpha=0.1, inherit.aes=FALSE
  ) +
  
  # Add true rate curve (TrueRate)
  geom_line(aes(y=TrueRate), color="darkgreen", size=1.2) +
  geom_point(aes(y=TrueRate), color="darkgreen", size=2) +
  
  labs(
    title="Simple Cluster Survival Test with Covariates",
    subtitle=paste0("Detected cluster: [", cs, ", ", ce, "]"),
    x="Time Index",
    y="Mortality Rate"
  ) +
  theme_minimal()

# If not yet installed:
# install.packages(c("Rcpp","ggplot2"))

library(Rcpp)
library(ggplot2)

# Step 1: Compile and load your current program (with runClusterSurvivalCovariates)
# sourceCpp("cluster_survival_covariates.cpp")

set.seed(20230908)

# 1. Generate more "natural-like" random data
N <- 30

# (a) At-risk population pop randomly between [80,150]
pop_r <- round(runif(N, min=80, max=150))

# (b) Multidimensional covariates (here, 3 examples):
#     Age (25~70), BMI (18~30), Gender (0/1)
Age    <- round(runif(N, min=25, max=70))
BMI    <- round(runif(N, min=18, max=30), 1)
Gender <- rbinom(N, 1, 0.45)
X      <- cbind(Age, BMI, Gender)    # 3 covariates

# (c) Set a range (e.g., 10~18) with higher mortality rates, lower elsewhere, plus random noise
base_p <- runif(N, 0.02, 0.05)   # Baseline mortality rate in [0.02,0.05]
for(i in 10:18){
  base_p[i] <- base_p[i] + 0.10 # Add 0.10 to form a hotspot
}
# Add some random fluctuation:
true_p <- pmax(pmin(base_p + rnorm(N, 0, 0.01), 0.99), 0.0001)

# (d) Generate number of deaths
deaths_r <- rbinom(N, size=pop_r, prob=true_p)

# (e) Organize into data.frame for viewing
df <- data.frame(
  Time      = 1:N,
  Pop       = pop_r,
  Deaths    = deaths_r,
  Observed  = deaths_r / pop_r,
  TrueRate  = true_p,
  Age       = Age,
  BMI       = BMI,
  Gender    = Gender
)
print(df)

# 2. Call the model
res <- runClusterSurvivalCovariates(deaths_r, pop_r, X, criterion_type=0) 
# 0 => BIC, 1 => AIC
print(res)

cs <- res$cs
ce <- res$ce
cat("Detected cluster: [", cs, ",", ce, "]\n")

# 3. Visualization
# Let the red rectangle only cover the data range (ymin, ymax), not the entire plot
ymin_val <- min(df$Observed, df$TrueRate)*0.9
ymax_val <- max(df$Observed, df$TrueRate)*1.1

ggplot(df, aes(x=Time)) +
  # (a) Bar chart showing Observed mortality rates
  geom_col(aes(y=Observed), fill="skyblue", alpha=0.6, width=0.8) +
  
  # (b) Highlight cluster range within y=[ymin_val, ymax_val]
  geom_rect(
    aes(
      xmin = cs - 0.4,
      xmax = ce + 0.4,
      ymin = ymin_val,
      ymax = ymax_val
    ),
    fill="red", alpha=0.2,
    inherit.aes=FALSE
  ) +
  
  # (c) Plot true rates curve (for demonstration or debugging)
  geom_line(aes(y=TrueRate), color="darkgreen", size=1) +
  geom_point(aes(y=TrueRate), color="darkgreen", size=2) +
  
  labs(
    title="Natural-like Test with Covariates (More Randomness)",
    subtitle=sprintf("Detected cluster: [%d, %d]", cs, ce),
    x="Time Index",
    y="Mortality Rate"
  ) +
  # Limit coordinate range
  coord_cartesian(ylim=c(0, max(df$Observed, df$TrueRate)*1.2)) +
  theme_minimal()
