## General Parameters
set.seed(12345677)

## birth and death process:
RUNS <- 10000
T_max <- 18.340862
alpha <- 0.018 # individual death rate
beta <- 0 # individual reproduction rate

for (run in seq(RUNS)){
  i = 1
  t = c(0)
  y = c(272)
  while(t[i] < T_max && y[i] > 0){
    i = i+1
    lambda = beta*y[i-1]
    mu = alpha*y[i-1]
    t[i] = t[i-1] + rexp(n=1, rate = mu+lambda) #sojourn time is exponentially distributed
    y[i] = y[i-1] + rbinom(n=1, prob = lambda/(lambda + mu), size = 1) *2 -1
  }
  if(run == 1){
    plot(t, y, t='s', 
         col = run, 
         main = "Pure Death Process Simulation Trend", 
         xlim = c(0, T_max), 
         ylim = c(0, 300),
         ylab = "Population Size",
         xlab = "Years of Survival",
         xaxt="n")
    axis(1, at = seq(0,18, by=2))
  } else {
    points(t, y, t='s', col = run)
  }
  #print(y[i])
}


for (run in seq(RUNS)){
  i = 1
  t = c(0)
  y = c(272)
  while(t[i] < T_max && y[i] > 0){
    i = i+1
    lambda = beta*y[i-1]
    mu = alpha*y[i-1]
    t[i] = t[i-1] + rexp(n=1, rate = mu+lambda) #sojourn time is exponentially distributed
    y[i] = y[i-1] + rbinom(n=1, prob = lambda/(lambda + mu), size = 1) *2 -1
  }
  if(run == 1){
    plot(t, y, t='s', 
         col = run, 
         main = "Pure Death Process: 15-18.34 Years", 
         xlim = c(15,18), 
         ylim = c(0, 300),
         ylab = "Population Size",
         xlab = "Years of Survival",
         xaxt="n")
    axis(1, at = seq(15,18, by=1))
  } else {
    points(t, y, t='s', col = run)
  }
  #print(y[i])
}


# par(mfrow=c(2,2))

# 0-5 years
axis(1, at = seq(0,5, by=1))
# 5-10 years
axis(1, at = seq(5,10, by=1))
# 10-15 years
axis(1, at = seq(10,15, by=1))
# 15-18 years
axis(1, at = seq(15,18, by=1))
# 0-18 years
axis(1, at = seq(0,18, by=2))


NKI_cleaned <- read.csv("~/UTH Stochastic Processes in Biostatistics I/PHD1950 Project/PHD1950 Data/NKI_cleaned.csv")

## birth and death process:
RUNS <- 10000
T_max <- NKI_cleaned$survival[272]
alpha <- 0.018 # individual death rate
beta <- 0 # individual reproduction rate

for (run in seq(RUNS)){
  i = 1
  t = c(0)
  y = c(272)
  while(t[i] < T_max && y[i] > 0){
    i = i+1
    lambda = beta*y[i-1]
    mu = alpha*y[i-1]
    t[i] = t[i-1] + rexp(n=1, rate = mu+lambda) #sojourn time is exponentially distributed
    y[i] = y[i-1] + rbinom(n=1, prob = lambda/(lambda + mu), size = 1) *2 -1
  }
  if(run == 1){
    plot(t, y, t='s', 
         col = run, 
         main = "Stochastic No Birth and Death Process", 
         xlim = c(0, T_max+1), 
         ylim = c(0, 300),
         ylab = "Population Size",
         xlab = "Years of Survival")
  } else {
    points(t, y, t='s', col = run)
  }
  #print(y[i])
}

# Not get the data to print as a vector. I tired many approaches for hours - limited by work efficiency
# This gets me the coverage interval across RUNS but not coverage probability
# hand to do coverage probability in excel
# Got it in R but coverage probability was always 100% so used Excel

# Code run in R to get the Kaplan Meier end of year populations

set.seed(12345677)

## birth and death process:
RUNS <- 100
T_max <- 19 # max time
alpha <- 0.018 # individual death rate
beta <- 0 # individual reproduction rate

for (run in seq(RUNS)){
  i = 1
  t = c(0)
  y = c(272)
  while(t[i] < T_max && y[i] > 0){
    i = i+1
    lambda = beta*y[i-1]
    mu = alpha*y[i-1]
    t[i] = t[i-1] + rexp(n=1, rate = mu+lambda) #sojourn time is exponentially distributed
    y[i] = y[i-1] + rbinom(n=1, prob = lambda/(lambda + mu), size = 1) *2 -1
  }
  if(run == 1){
    plot(t, y, t='s', 
         col = run, 
         main = "Stochastic No Birth and Death Process", 
         xlim = c(0, T_max+1), 
         ylim = c(0, 300),
         ylab = "Population Size",
         xlab = "Years of Survival")
  } else {
    points(t, y, t='s', col = run)
  }
  print(y[i])
}




## General Parameters
set.seed(12345677)

## birth and death process:
RUNS <- 100
T_max <- 18.340862*365+4 # max time
alpha <- 0.018/365 # individual death rate
beta <- 0 # individual reproduction rate

for (run in seq(RUNS)){
  i = 1
  t = c(0)
  y = c(272)
  while(t[i] < T_max && y[i] > 0){
    i = i+1
    lambda = beta*y[i-1]
    mu = alpha*y[i-1]
    t[i] = t[i-1] + rexp(n=1, rate = mu+lambda) #sojourn time is exponentially distributed
    y[i] = y[i-1] + rbinom(n=1, prob = lambda/(lambda + mu), size = 1) *2 -1
  }
  if(run == 1){
    plot(t, y, t='s', 
         col = run, 
         main = "Stochastic No Birth and Death Process", 
         xlim = c(0, T_max+1), 
         ylim = c(0, 300),
         ylab = "Population Size",
         xlab = "Days of Survival")
  } else {
    points(t, y, t='s', col = run)
  }
  #print(y[i])
}




# can not be used with real data

# Code I am trying above but never got to work
pop_run <- print(y[i]) 
mean_pop <- y[i]
# get the end result of the population over the # of runs
y_pop <- data.frame(mean_pop)
view(y_pop)


my_vector <- vector("numeric", RUNS)
for(i in 1:RUNS) my_vector[y[i]] <- i
my_vector

ci <- function(i){
  lo <- m[i] - 1.96 * (m[i]*(1 - m[i])/n)^0.5
  hi <- m[i] + 1.96 * (m[i]*(1 - m[i])/n)^0.5
  c(lo, hi)
}

p1 <- 0.01
p2 <- 0.05
p3 <- 0.1
p4 <- 0.3
p5 <- 0.5
p6 <- 0.7
p7 <- 0.9
p8 <- 0.95
p9 <- 0.99
p10 <- 0.9999
p <- c(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)

m <- length(p)
n <- 1000
Runs <- 10000
CI <- array(NA, dim = c(m, 3, Runs))
set.seed (123467)    # do this just once outside the loop
for (i in 1:Runs){
  x <- matrix(rbinom(10*n, 1, rep(p, each = n)), ncol = 10)
  m <- colMeans(x)
  CI[, 1:2, i] <- t(sapply(seq_along(m), ci))
  CI[, 3, i] <- CI[, 2, i] - CI[, 1, i]
}
dimnames(CI) <- list(c("p=0,01","p=0.05","p=0.1","p=0,3","p=0,5","p=0,7", "p=0,9","p=0,95","p=0,99", "p=0,9999"),
                     c("lower.bound","upper.bound", "width"),
                     sprintf("R%05d", 1:Runs))

CI[,, 200]  # example result of loop

bad <- array(NA, dim = c(length(p), Runs))
for(i in 1:Runs){
  bad[, i] <- as.integer(!(CI[, 1, i] < p & p < CI[, 2, i]))
}

prop_bad <- rowMeans(bad)
names(prop_bad) <- c("p=0,01","p=0.05","p=0.1","p=0,3","p=0,5","p=0,7", "p=0,9","p=0,95","p=0,99", "p=0,9999")
prop_bad
