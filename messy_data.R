library(tidyverse)
library(nlme)

set.seed(123)

n_people <- 60
n_visits <- 6

dat <- data.frame(
  id = factor(rep(seq_len(n_people), each = n_visits)),
  week = rep(0:5, times = n_people),
  product = factor(
    rep(rep(c("A", "B"), each = 30), each = n_visits)
  )
)

#participant deviations from the average starting level
person_intercept <- rep(
  rnorm(n_people, mean = 0, sd = 8),
  each = n_visits
)

#participant deviations from their product's average slope
person_slope <- rep(
  rnorm(n_people, mean = 0, sd = 2),
  each = n_visits
)

#nearby measurements within a participant are more correlated
#arima.sim generates n_visits correlated errors for one participant
#replicate repeats this for n_people
#e(t) = 0.5*e(t-1) + u(t), where u(t)~N(0,sd=3)
measurement_error <- as.vector(
  replicate(
    n_people,
    arima.sim(
      model = list(ar = 0.5),
      n = n_visits,
      sd = 3
    )
  )
)

#both products improve over time
#with product B improves an additional 1.5 points per week
dat$satisfaction <- with(
  dat,
  50 +
    person_intercept +
    (1 + 1.5 * (product == "B") + person_slope) * week +
    measurement_error
)

head(dat)

set.seed(456)

dat_original <- dat %>%
  select(id, week, product, satisfaction)

messy_dat <- dat_original %>%
  mutate(
    id = paste0(" P", as.character(id), " "),
    
    week = paste("Week", week),

#replace "A" with either "a" or " Product A "; same for "B"    
    product = case_when(
      product == "A" ~ sample(
        c("A", "a", " Product A "),
        n(), replace = TRUE
      ),
      product == "B" ~ sample(
        c("B", "b", " Product B "),
        n(), replace = TRUE
      )
    ),
    
    satisfaction = paste0(
      round(satisfaction, 2),
      " points"
    )
  )

#introduce missing scores in 12 different rows
missing_rows <- sample(seq_len(nrow(messy_dat)), 12)

messy_dat$satisfaction[missing_rows] <- NA_character_

#introduce two additional problems in different rows
remaining_rows <- setdiff(
  seq_len(nrow(messy_dat)),
  missing_rows
)

problem_rows <- sample(remaining_rows, 2)

#known data entry error
messy_dat$satisfaction[problem_rows[1]] <- "999 points"

#missing visit information
messy_dat$week[problem_rows[2]] <- "unknown"

#add five exact duplicate records
messy_dat <- bind_rows(
  messy_dat,
  messy_dat %>% slice_sample(n = 5)
)

glimpse(messy_dat)
head(messy_dat)

####################
#inspect the problems
messy_dat %>% count(product)

#visit labels
messy_dat %>% count(week)

#missing values in each column
messy_dat %>%
  summarise(across(everything(), ~ sum(is.na(.x))))

#number of exact duplicate rows
sum(duplicated(messy_dat))

summary(messy_dat)


#####################
#clean the data
clean_dat <- messy_dat %>%
  distinct() %>%                     #remove exact duplicates
  mutate(
    id = str_trim(id),
    
    product = str_to_upper(str_trim(product)),
    product = str_remove(product, "^PRODUCT\\s+"),
    
    # Keep the original text for checking
    week_raw = week,
    satisfaction_raw = satisfaction,
    
    week = parse_number(
      week,
      na = c("", "NA", "unknown")
    ),
    
    satisfaction = parse_number(satisfaction),
    
    invalid_score = !is.na(satisfaction) &
      satisfaction == satisfaction[which(!is.na(satisfaction) & satisfaction > 100)],
    
    # Treat the known entry error as missing
    satisfaction = replace(
      satisfaction,
      invalid_score,
      NA_real_
    ),
    
    id = factor(id),
    product = factor(product, levels = c("A", "B"))
  ) %>%
  arrange(id, week)

#invalid records
flagged_records <- clean_dat %>%
  filter(
    invalid_score |
      is.na(satisfaction) |
      is.na(week) |
      is.na(product)
  )

flagged_records %>%
  select(
    id, product,
    week_raw, satisfaction_raw,
    invalid_score
  )

#check the longitudinal structure, where 
#each participant should have at most 1 one record/week
#and product assignment should remain constant
#should return zero rows
duplicate_visits <- clean_dat %>%
  filter(!is.na(week)) %>%
  count(id, week) %>%
  filter(n > 1)

duplicate_visits

#should return zero rows
product_changes <- clean_dat %>%
  group_by(id) %>%
  summarise(
    n_products = n_distinct(product, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_products > 1)

product_changes

#stop if either problem remains
stopifnot(
  nrow(duplicate_visits) == 0,
  nrow(product_changes) == 0
)

#exclude unusable visits, while retaining the participant's valid visits
analysis_dat <- clean_dat %>%
  filter(
    !is.na(id),
    !is.na(product),
    !is.na(week),
    !is.na(satisfaction)
  ) %>%
  arrange(id, week)

#report how many records remain
tibble(
  stage = c(
    "Original data",
    "Messy data",
    "After removing duplicates",
    "Available for analysis"
  ),
  rows = c(
    nrow(dat_original),
    nrow(messy_dat),
    nrow(clean_dat),
    nrow(analysis_dat)
  )
)

#number of usable visits per participant
analysis_dat %>%
  count(id, name = "usable_visits") %>%
  count(usable_visits, name = "participants")


#fit the models
#random intercept and slope, independent residual errors
model_no_ar <- lme(
  fixed = satisfaction ~ product * week,
  random = ~ week | id,
  data = analysis_dat,
  method = "REML",
  na.action = na.fail
)

#add residual AR(1) correlation
model_ar <- update(
  model_no_ar,
  correlation = corAR1(form = ~ week | id)
)

#compare covariance structures using the same records
anova(model_no_ar, model_ar)

summary(model_ar)
VarCorr(model_ar)

#estimate each product’s average weekly slope and compare them
library(emmeans)

product_slopes <- emtrends(
  model_ar,
  ~ product,
  var = "week"
)

#slopes and 95% confidence intervals
confint(product_slopes)

#Product B slope minus Product A slope
summary(
  contrast(product_slopes, method = "revpairwise"),
  infer = c(TRUE, TRUE)
)

#plot the cleaned observations and fitted lines
analysis_dat$fitted <- predict(
  model_ar,
  newdata = analysis_dat,
  level = 1
)

ggplot(
  analysis_dat,
  aes(x = week, color = product, group = id)
) +
  geom_point(
    aes(y = satisfaction),
    alpha = 0.3,
    size = 1
  ) +
  geom_line(
    aes(y = fitted),
    alpha = 0.6
  ) +
  facet_wrap(~ product) +
  scale_x_continuous(breaks = 0:5) +
  labs(
    title = "Cleaned data and participant fitted lines",
    x = "Week",
    y = "Satisfaction",
    color = "Product"
  ) +
  theme_minimal()
