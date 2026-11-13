
library(MASS)
library(nlme)

simulate_crossover <- function(n,
                               delta = 3,
                               sigma = 8,
                               rho = 0.5) {
  
  stopifnot(n %% 2 == 0)
  
  #balanced AB and BA sequences
  sequence <- rep(c("AB", "BA"), each = n / 2)
  
  dat <- expand.grid(
    day = 1:3,
    period = 1:2,
    id = seq_len(n)
  )
  
  dat <- dat[order(dat$id, dat$period, dat$day), ]
  
  dat$sequence <- factor(sequence[dat$id])
  
  dat$treatment <- ifelse(
    (dat$sequence == "AB" & dat$period == 1) |
      (dat$sequence == "BA" & dat$period == 2),
    "A", "B"
  )
  
  dat$treatment <- factor(dat$treatment)
  
  dat$visit <- (dat$period - 1) * 3 + dat$day
  
  # Correlation across all six observations
  R <- outer(1:6, 1:6,
             function(x, y) rho^abs(x - y))
  
  Sigma <- sigma^2 * R
  
  errors <- MASS::mvrnorm(
    n = n,
    mu = rep(0, 6),
    Sigma = Sigma
  )
  
  # Expected treatment difference = delta
  dat$y <- 40 +
    delta * (dat$treatment == "B") +
    0.5 * (dat$period == 2) +
    0.2 * dat$day +
    as.vector(t(errors))
  
  dat
}

estimate_power_crossover <- function(n,
                                     nsim = 1000) {
  
  pvalues <- replicate(nsim, {
    
    dat <- simulate_crossover(n)
    
    tryCatch({
      
      full <- gls(
        y ~ treatment + factor(period) +
          sequence + day,
        correlation = corAR1(form = ~ visit | id),
        data = dat,
        method = "ML"
      )
      
      reduced <- gls(
        y ~ factor(period) + sequence + day,
        correlation = corAR1(form = ~ visit | id),
        data = dat,
        method = "ML"
      )
      
      anova(reduced, full)$`p-value`[2]
      
    }, error = function(e) NA_real_)
    
  })
  
  c(
    power = mean(pvalues < 0.05, na.rm = TRUE),
    failures = sum(is.na(pvalues))
  )
}

set.seed(123)

sample_sizes <- c(20, 30, 40, 60, 80)

results_crossover <- t(sapply(
  sample_sizes,
  estimate_power_crossover,
  nsim = 1000
))

data.frame(
  total_n = sample_sizes,
  results_crossover
)
