#calculate power using simulation in R
#treatment vs control longitudinal study


library(MASS)
library(nlme)

simulate_parallel <- function(n_per_group,
                              delta = 3,    #effect we want to detect at the end of the study
                              sigma = 8,    #standard deviation at each timepoint
                              rho = 0.5) {  #correlation between adjacent observations
  
  times <- 0:4
  n <- 2 * n_per_group
  
  group <- rep(c(0, 1), each = n_per_group)
  
  # AR(1) covariance matrix
  R <- outer(times, times,
             function(x, y) rho^abs(x - y))
  
  Sigma <- sigma^2 * R
  
  #generate correlated residuals for each subject
  errors <- MASS::mvrnorm(    #mvrnorm generates random obs from a multivariable normal distribution, which can be correlated
    n = n,
    mu = rep(0, length(times)),
    Sigma = Sigma
  )
  
  dat <- expand.grid(
    time = times,
    id = seq_len(n)
  )
  
  dat <- dat[order(dat$id, dat$time), ]
  
  dat$group <- group[dat$id]
  
  #difference at final visit = delta
  dat$y <- 40 +
    delta * dat$group * dat$time / 4 +
    as.vector(t(errors))
  
  dat
}

estimate_power_parallel <- function(n_per_group,
                                    nsim = 1000) {
  
  pvalues <- replicate(nsim, {
    
    dat <- simulate_parallel(n_per_group)
    
    tryCatch({
      
      full <- gls(
        y ~ group * time,
        correlation = corAR1(form = ~ time | id),
        data = dat,
        method = "ML"
      )
      
      reduced <- gls(
        y ~ group + time,
        correlation = corAR1(form = ~ time | id),
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

sample_sizes <- c(30, 50, 80, 120)

results_parallel <- t(sapply(
  sample_sizes,
  estimate_power_parallel,
  nsim = 1000
))

data.frame(
  n_per_group = sample_sizes,
  total_n = 2 * sample_sizes,
  results_parallel
)

#for evaluating sensitivity to the assumptions,
#repeat the analysis for several treatment differences (2,3,4), 
#lag 1 correlations (0.2,0.5,0.8),
#standard deviations (6,8,10),
#measurement visits (3,5,7).