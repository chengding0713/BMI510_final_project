# BMI510_final_project

Author
------------
Cheng Ding (cheng.ding2@emory.edu)

Functions
---------

bmi510package includes the following functions:

-   `rando()`: A wrapper around `sample()` for randomly sampling atomic vectors or dataframe-like objects.
-   `is_min()` and `is_max()`: Functions to identify the minimum or maximum values in an atomic vector.
-   `rep_mat()`: A port of the `repmat.m` function from MATLAB, used for replicating matrix rows or columns.
-   `classes()`: Returns a character vector containing the classes of each variable in a tibble.
-   `df_scale()`: Scales the numeric variables in a tibble with optional centering and scaling.
-   `log_likelihood_*()`: A set of functions to calculate log-likelihoods under various distributions (normal, uniform, chi-squared, F, and t).
-   `sensitivity()`, `specificity()`, `precision()`, `recall()`, `accuracy()`, and `f1()`: Functions to calculate various performance metrics for binary classifiers.
-   `minimum_n_per_group()`: Returns the minimum sample size per group needed for a two-sample t-test, based on the expected Cohen's d and desired statistical power.
-   `r2()`: Calculates the R-squared statistic between predicted and ground truth continuous variables.
-   `adj_R2()`: Calculates the adjusted R-squared statistic between predicted and ground truth continuous variables, accounting for the number of model parameters.

Examples
--------
<pre>
```R
`# Randomly sample rows from a data.frame
data(mtcars)
sampled_rows <- rando(mtcars, n = 5, replace = FALSE)

# Calculate log-likelihood under normal distribution
x <- rnorm(100, mean = 0, sd = 1)
log_likelihood_norm(x, mean = 0, sd = 1)

# Evaluate classifier performance
pred <- factor(c(1, 0, 1, 1, 0))
truth <- factor(c(1, 0, 1, 0, 0))
accuracy(pred, truth)`
```
</pre>
