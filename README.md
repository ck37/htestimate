# htestimate

Htestimate calculates unbiased estimates of treatment effects from randomized trials when the random assignment is correlated across units, using the Horvitz-Thompson estimator (Särndal et al. 2003, section 2.8). Standard approaches to RCT evaluation (difference in means and regression) are generally biased under clustered randomization (Middleton 2008) or under rerandomization (Morgan & Rubin 2012), for example. In addition to the treatment effect the package produces a standard error and p-value of that effect estimate. Differences in outcome totals rather than means can also be produced. Any number of experimental arms/conditions are allowed.

This package is currently under active development so bug reports and feature requests are encouraged.

## Install

Install directly from github using devtools (install.packages("devtools") if you don't already have devtools):
```{r}
library(devtools)
install_github("ck37/htestimate")
library(htestimate)
```

## Requirements

R packages: dplyr

## Examples

```{r}
# Example using data from RI package.
y <- c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
Z <- c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
# Generate 10,000 random permutations of the assignment vector.
perms = ri::genperms(Z, maxiter=10000)
# Estimate the probability of assignment for each unit and assignment level.
prob_matrix = createProbMatrix(perms)
# Estimate the treatment effect using Horvitz-Thompson.
htestimate(y, Z, contrasts = c(-1, 1), prob_matrix = prob_matrix)
```

## References

Aronow, P. M., & Middleton, J. A. (2013). A class of unbiased estimators of the average treatment effect in randomized experiments. Journal of Causal Inference, 1(1), 135-154.

Middleton, J. A. (2008). Bias of the regression estimator for experiments using clustered random assignment. Statistics & Probability Letters, 78(16), 2654-2659.

Morgan, K. L., & Rubin, D. B. (2012). Rerandomization to improve covariate balance in experiments. The Annals of Statistics, 40(2), 1263-1282.

Särndal, C. E., Swensson, B., & Wretman, J. (2003). Model assisted survey sampling. Springer Science & Business Media.
