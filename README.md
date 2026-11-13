# Longitudinal Satisfaction Analysis in R

This project explores repeated satisfaction measurements using simulated data, deliberate data quality problems, and mixed effects models. It demonstrates how to clean longitudinal records, model participant differences, assess residual correlation, and communicate effects with ggplot2.

All observations are synthetic. The analysis is an educational example rather than evidence about actual products or consumers.

## Research question

Does satisfaction change over time, and does the average rate of change differ between Product A and Product B?

## Simulated study

The study contains 60 participants measured weekly from week 0 through week 5, giving 360 observations before introducing data quality problems. Thirty participants receive each product. Product labels are assigned by participant number in the simulation rather than randomized.

| Variable | Description |
| --- | --- |
| `id` | Participant identifier |
| `week` | Time since the first visit, from 0 through 5 |
| `product` | Product A or Product B |
| `satisfaction` | Simulated continuous satisfaction score |

The final simulation includes participant variation in both starting satisfaction and weekly slopes. Earlier examples used random intercepts alone.

The generating model is:

$$
Y_{it} = 50 + b_{0i} + [1 + 1.5B_i + b_{1i}]t + e_{it}
$$

Here, $B_i$ equals 1 for Product B and 0 for Product A. Participant intercept deviations have standard deviation 8, and slope deviations have standard deviation 2. These deviations are generated independently.

Product A has an average generating slope of 1 point per week. Product B has an average generating slope of 2.5 points per week. The generating difference in average slopes is therefore 1.5 points per week. These are simulation parameters, not fitted results.

Errors follow an AR(1) process with correlation parameter 0.5. The innovation standard deviation is 3, corresponding to a stationary error standard deviation of approximately 3.46.

## Data quality exercise

A copy of the simulated data was deliberately modified to include:

* Extra whitespace in participant identifiers.
* Inconsistent product labels such as `A`, `a`, and ` Product A `.
* Visit labels stored as text, such as `Week 3`.
* Satisfaction stored as text, such as `52.4 points`.
* Twelve missing satisfaction scores.
* One known entry error recorded as `999 points`.
* One visit recorded as `unknown`.
* Five exact duplicate records.

The original data are retained as `dat_original`. The modified data are stored as `messy_dat`.

## Cleaning workflow

1. Remove exact duplicate records with `distinct()`.
2. Trim surrounding whitespace with `str_trim()`.
3. Standardize product labels using `str_to_upper()` and `str_remove()`.
4. Extract numeric weeks and satisfaction scores with `parse_number()`.
5. Preserve the original week and satisfaction text for review.
6. Flag the known entry error and replace its score with a missing value.
7. Check for conflicting records within participant and week.
8. Check that product assignment remains constant within participant.
9. Create `analysis_dat` using visits with complete model variables.

Exact duplicates can be removed in this exercise because their origin is known. Conflicting visit records would require investigation. The value 999 is treated as an error because it was deliberately inserted, not because all unusual observations should be discarded.

Participants retain their other valid visits when a measurement is missing. Original week values are preserved so that missing visits do not change elapsed time.

## R packages

```r
install.packages(c("tidyverse", "nlme", "emmeans"))

library(tidyverse)
library(nlme)
library(emmeans)
```

## Mixed effects model

The main specification includes product, week, and their interaction, participant random intercepts and slopes, and an AR(1) residual structure.

```r
model_ar <- lme(
  fixed = satisfaction ~ product * week,
  random = ~ week | id,
  correlation = corAR1(form = ~ week | id),
  data = analysis_dat,
  method = "REML",
  na.action = na.fail
)
```

| Component | Interpretation |
| --- | --- |
| Intercept | Average starting satisfaction for Product A |
| `productB` | Difference between products at week 0 |
| `week` | Average weekly change for Product A |
| `productB:week` | Difference in average weekly slopes between products |
| Random intercept | Participant deviation in starting satisfaction |
| Random slope | Participant deviation in weekly change |
| AR(1) | Additional residual dependence within participant |

The random intercept and slope specification also estimates their correlation. A random slope does not automatically require AR(1). These represent different sources of dependence.

## Assessing residual correlation

A model with random intercepts and slopes but independent errors was compared with the AR(1) model:

```r
model_no_ar <- update(model_ar, correlation = NULL)

anova(model_no_ar, model_ar)

plot(
  ACF(model_ar, resType = "normalized"),
  alpha = 0.05
)
```

The comparison uses identical fixed effects and observations, so REML likelihood comparisons are appropriate. AIC, the likelihood ratio test, and residual diagnostics should be interpreted together.

Residual ACF plots need care with only six visits. Estimating participant intercepts and slopes can introduce negative residual correlations. Approximate plotted limits do not fully account for that estimation. No final numerical conclusion about adding AR(1) is documented here because its comparison output was not supplied.

## Assessing curvature

A quadratic time model was compared with the linear time model using ML because their fixed effects differ:

```r
model_linear <- update(model_ar, method = "ML")

model_quadratic <- update(
  model_linear,
  fixed = satisfaction ~ product * (week + I(week^2))
)

anova(model_linear, model_quadratic)
```

The analysis produced the following supplied comparison:

| Model | Parameters | AIC | BIC | Log likelihood |
| --- | --- | --- | --- | --- |
| Linear | 9 | 2064.709 | 2099.327 | −1023.355 |
| Quadratic | 11 | 2068.416 | 2110.727 | −1023.208 |

The likelihood ratio statistic was 0.2930, with 2 additional parameters and a p value of 0.8637.

The quadratic terms provided little improvement. Both AIC and BIC favored the linear model. This supports retaining linear time for this example and is consistent with the generating process.

These values describe the run shared during the analysis. Regenerating data, changing exclusions, or changing model specifications can change the results.

A quadratic trajectory is nonlinear in time but remains linear in its coefficients. A true nonlinear mixed effects model, such as exponential growth toward a plateau, was discussed but was not fitted.

## Final estimation and interpretation

Refit the selected fixed effects structure using REML:

```r
model_final <- update(model_linear, method = "REML")

summary(model_final)
VarCorr(model_final)
```

Estimated product means at each week:

```r
estimated_means <- emmeans(
  model_final,
  ~ product | week,
  at = list(week = 0:5)
)

confint(estimated_means)
```

Product B minus Product A at each week:

```r
product_effects <- contrast(
  estimated_means,
  method = "revpairwise",
  adjust = "none"
)

confint(product_effects)
```

Average weekly slopes and their difference:

```r
product_slopes <- emtrends(
  model_final,
  ~ product,
  var = "week"
)

confint(product_slopes)

summary(
  contrast(product_slopes, method = "revpairwise"),
  infer = c(TRUE, TRUE)
)
```

A product difference at a particular week is measured in satisfaction points. A slope difference is measured in satisfaction points per week. These answer different questions.

## Visualizations and diagnostics

The analysis developed ggplot2 code for:

* Observed average satisfaction by product and week.
* Model estimated mean trajectories with 95% confidence intervals.
* Product differences at each week, labeled with estimated effect sizes.
* Participant fitted lines including random effects.
* Participant fitted curves from the exploratory quadratic model.

Participant predictions use `predict(model_final, newdata = analysis_dat, level = 1)`. Population predictions use `level = 0`.

Error bars for estimated means describe uncertainty in population means, not the spread of individual observations. Product contrasts use pointwise intervals without adjustment across weeks.

Diagnostics discussed include residuals versus fitted values, residual normal quantile plots, and the ACF of normalized residuals.

## Reproducibility and limitations

The simulation uses seed 123, and the deliberate data corruption uses seed 456. Run the simulation, corruption, cleaning, and modeling code in that order. Save `sessionInfo()` with analysis outputs to record the R environment.

The README documents the analysis developed in this project. Its code excerpts assume the corresponding data preparation steps have already created `analysis_dat`.

The simulated outcome is continuous and unbounded. It does not enforce the limits of a real satisfaction scale. Product assignment and missingness mechanisms are simplified for method illustration.

Mixed models can use incomplete participant records without filling in missing outcomes. Valid inference still depends on appropriate model assumptions and an ignorable missingness mechanism, typically missing at random conditional on the modeled and observed information.

The primary demonstrated model selection result concerns linear versus quadratic time. Fitted treatment effect estimates should be reported from the actual model output rather than substituted with the known simulation parameters.

