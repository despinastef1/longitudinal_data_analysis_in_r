#example using a linear mixed effects model for analyzing 
#repeated measurements from the same participants

#suppose 60 participants are assigned to Product A or Product B
#the satisfaction scores are measured weekly for 6 weeks

#Question: Does satisfaction change differently over time for the two products?

#Step 1: create example data
library(nlme)
library(emmeans)
library(ggplot2)

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

#plot average trajectories
average_scores <- aggregate(
  satisfaction ~ week + product,
  data = dat,
  FUN = mean
)
average_scores <- as.data.frame(average_scores)

ggplot(data=average_scores, aes(x = week, y=satisfaction,color=product)) +
  geom_point() +
  geom_line() +
  labs(
    x = "Week",
    y = "Mean satisfaction",
    color = "Product"
  ) +
  theme_minimal()

#fit the longitudinal model
model <- lme(
  fixed = satisfaction ~ product * week,  #to check for changes in scores that differ over time
  random = ~ week| id,  #respondents are random
  correlation = corAR1(form = ~ week | id),  #this accounts for repeated measurements
  data = dat,
  method = "REML"
)

summary(model)

#average effects
fixef(model)

#tests for the fixed effects
anova(model, type = "marginal")

#confidence intervals for the fixed effects
intervals(model, which = "fixed")

#participant deviations from the average intercept and slope
ranef(model)

#estimated variability in participant intercepts and slopes
VarCorr(model)

#estimated residual correlation between adjacent visits
coef(
  model$modelStruct$corStruct,
  unconstrained = FALSE
)

#estimate mean satisfaction for each product at each week
estimated_means <- emmeans(
  model,
  ~ product | week,
  at = list(week = 0:5)
)

means_df <- as.data.frame(confint(estimated_means))
means_df

#estimated trajectories with 95% confidence intervals for mean satisfaction
ggplot(
  means_df,
  aes(
    x = week,
    y = emmean,
    color = product,
    group = product
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.12
  ) +
  scale_x_continuous(breaks = 0:5) +
  labs(
    title = "Estimated satisfaction over time",
    subtitle = "Error bars show 95% confidence intervals",
    x = "Week",
    y = "Estimated mean satisfaction",
    color = "Product"
  ) +
  theme_minimal(base_size = 13)

#effect sizes
#with product levels A, B, this calculates B minus A
product_effects <- contrast(
  estimated_means,
  method = "revpairwise",
  adjust = "none"
)

effects_df <- as.data.frame(confint(product_effects))
effects_df

ggplot(effects_df, aes(x = week, y = estimate)) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    color = "gray50"
  ) +
  geom_line(color = "purple", linewidth = 1) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.12,
    color = "purple"
  ) +
  geom_point(size = 3, color = "purple") +
  geom_text(
    aes(
      y = upper.CL,
      label = sprintf("%.2f", estimate)
    ),
    vjust = -0.7,
    show.legend = FALSE
  ) +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(
    expand = expansion(mult = c(0.1, 0.2))
  ) +
  labs(
    title = "Product effect at each week",
    subtitle = "Labels show estimated differences in satisfaction points",
    x = "Week",
    y = "Product B minus Product A",
    caption = "Error bars show pointwise 95% confidence intervals."
  ) +
  theme_minimal(base_size = 13)

#plot fitted line for each participant
#use predict(model, level = 1) to obtain fitted values that include each participant’s random intercept
#group = id draws a separate line for each participant

dat$fitted <- predict(model, newdata = dat, level = 1)

ggplot(
  dat,
  aes(x = week, color = product, group = id)
) +
  geom_point(
    aes(y = satisfaction),
    alpha = 0.3,
    size = 1
  ) +
  geom_line(
    aes(y = fitted),
    alpha = 0.6,
    linewidth = 0.7
  ) +
  facet_wrap(~ product) +
  labs(
    title = "Fitted satisfaction for each participant",
    x = "Week",
    y = "Satisfaction",
    color = "Product"
  ) +
  theme_minimal(base_size = 13)

#check model assumptions
#look for changing residual spread or systematic patterns
plot(model)

#check approximate normality of residuals
qqnorm(resid(model, type = "normalized"))
qqline(resid(model, type = "normalized"))

#check remaining correlation within participants (after accounting for the fitted model)
plot(ACF(model, resType = "normalized"))

#now, fit the model without the AR(1) structure
#this helps to check if AR(1) structure is supported
dat <- dat[order(dat$id, dat$week), ]

model_no_ar <- lme(
  fixed = satisfaction ~ product * week,
  random = ~ week | id,
  data = dat,
  method = "REML"
)

acf_before <- ACF(
  model_no_ar,
  resType = "normalized",
  maxLag = 5
)

print(acf_before)
plot(acf_before, alpha = 0.05)

#add AR(1) and compare the models
model_ar <- update(
  model_no_ar,
  correlation = corAR1(form = ~ week | id)
)

anova(model_no_ar, model_ar)

coef(
  model_ar$modelStruct$corStruct,
  unconstrained = FALSE
)
acf_after <- ACF(
  model_ar,
  resType = "normalized",
  maxLag = 5
)

plot(acf_after, alpha = 0.05)

#try a curved trajectory model
#straight line model
model_linear <- lme(
  fixed = satisfaction ~ product * week,
  random = ~ week | id,
  correlation = corAR1(form = ~ week | id),
  data = analysis_dat,
  method = "REML"
)

#curved trajectory model
model_quadratic <- update(
  model_linear,
  fixed = satisfaction ~ product * (week + I(week^2))
)

summary(model_quadratic)

#does adding curvature improve the fit?
anova(model_linear, model_quadratic)
