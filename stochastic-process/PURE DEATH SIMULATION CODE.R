library(MASS)

# Set simulation parameters
num_simulations <- 10000 # number of simulations to run
num_years <- 20 # number of years to simulate
initial_population <- 10000 # initial population size 
# population = 272 for us
true_death_rate <- 0.02 # true annual death rate
alpha <- 0.05 # significance level for confidence interval

# Generate sample data using Poisson distribution
deaths <- matrix(nrow=num_simulations, ncol=num_years) # create matrix to store deaths
for (i in 1:num_simulations) {
  deaths[i, 1] <- 0 # no deaths in the first year
  for (j in 2:num_years) {
    deaths[i, j] <- rpois(1, true_death_rate * initial_population) # simulate deaths using Poisson distribution
  }
}

# Estimate death rate using maximum likelihood estimation
est_death_rate <- apply(deaths, 1, function(x) {
  fit <- fitdistr(x, "Poisson")
  fit$estimate
})/10000

# Calculate coverage probability for a 95% confidence interval
lower_bound <- est_death_rate - qnorm(1 - alpha/2) * sqrt(est_death_rate / initial_population)
upper_bound <- est_death_rate + qnorm(1 - alpha/2) * sqrt(est_death_rate / initial_population)
coverage_prob <- mean((lower_bound <= true_death_rate) & (true_death_rate <= upper_bound))

# Print results
cat("Coverage probability:", coverage_prob)

# Summary Statistic
summary(deaths)
