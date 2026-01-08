set.seed(1234567) # cov prob = 98.84

library(MASS)

# Set simulation parameters
num_simulations <- 10000 # number of simulations to run
num_years <- 18.340862 # number of years to simulate
initial_population <- 272 # initial population size
true_death_rate <- 0.26 # true annual death rate
alpha <- 0.05 # significance level for confidence interval

# Generate sample data using Poisson distribution
deaths <- matrix(nrow=num_simulations, ncol=num_years) # create matrix to store deaths
for (i in 1:num_simulations) {
  deaths[i, 1] <- 0 # no deaths in the first year
  for (j in 2:num_years) {
    deaths[i, j] <- rpois(1, true_death_rate * initial_population) # simulate deaths using Poisson distribution
  }
}
# popdrop <- initial_population - (deaths[1:(i-1)]) 
# Estimate death rate using maximum likelihood estimation
est_death_rate <- apply(deaths, 1, function(x) {
  fit <- fitdistr(x, "Poisson")
  fit$estimate
})/(initial_population-deaths[i-1, j-1])

# Calculate coverage probability for a 95% confidence interval
lower_bound <- est_death_rate - qnorm(1 - alpha/2) * sqrt(est_death_rate / (initial_population-deaths[i-1, j-1]))
upper_bound <- est_death_rate + qnorm(1 - alpha/2) * sqrt(est_death_rate / (initial_population-deaths[i-1, j-1]))
coverage_prob <- mean((true_death_rate >= lower_bound) & (true_death_rate <= upper_bound))

# Print results
cat("Coverage probability:", coverage_prob*100)

# Notes:
# Cov prob = 95.77 at 5 Years
# Cov prob = 99.49 at 10 Years
# Cov prob = 99.95 at 15 Years
# Cov prob = 98.84 at 18.340862 Years
