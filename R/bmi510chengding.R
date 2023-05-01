library(pwr)


#' Randomly sample elements or rows from a vector or data frame
#'
#' This function is a wrapper around the \code{sample} function that can handle both atomic vectors and data frames as input. It randomly samples either \code{n} elements from an atomic vector or \code{n} rows from a data frame, with or without replacement.
#'
#' @param x A vector or data frame to sample from
#' @param n The number of samples or rows to select (default = 1)
#' @param replace Logical. Should sampling be done with replacement? (default = TRUE)
#'
#' @return A vector or data frame containing the sampled elements or rows
#'
#' @export
rando = function(x, n=1, replace=T){
  if (is.atomic(x)) {
    return(sample(x, n, replace=replace))
  } else if (is.data.frame(x)) {
    return(x[sample(nrow(x), n, replace=replace), ])
  }
}


#' This function accepts an atomic vector \code{x} and returns a logical vector with \code{TRUE} values where \code{x} equals its minimum value.
#'
#' @param x An atomic vector to check
#' @param na.rm Logical. Should missing values be removed? (default = TRUE)
#'
#' @return A logical vector with \code{TRUE} values where \code{x} equals its minimum value
#'
#' @export
is_min <- function(x, na.rm = TRUE) {
  min_value <- min(x, na.rm = na.rm)
  return(x == min_value)
}

#' Check if values in a vector are equal to the maximum value
#'
#' This function accepts an atomic vector \code{x} and returns a logical vector with \code{TRUE} values where \code{x} equals its maximum value.
#'
#' @param x An atomic vector to check
#' @param na.rm Logical. Should missing values be removed? (default = TRUE)
#'
#' @return A logical vector with \code{TRUE} values where \code{x} equals its maximum value
#'
#' @export
is_max <- function(x, na.rm = TRUE) {
  max_value <- max(x, na.rm = na.rm)
  return(x == max_value)
}

#' Replicate a matrix or data frame by repeating rows or columns
#'
#' This function replicates a matrix or data frame \code{x} by repeating its rows \code{M} times or its columns \code{N} times, depending on the parameters passed.
#'
#' @param x A matrix or data frame to replicate
#' @param M The number of times to repeat the rows of \code{x} (default = 1)
#' @param N The number of times to repeat the columns of \code{x} (default = 1)
#'
#' @return A matrix with rows and/or columns replicated M and/or N times
#'
#' @export
rep_mat <- function(x, M = 1, N = 1) {
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("Invalid input: x must be a matrix or dataframe.")
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  if (rows == 0 || cols == 0) {
    stop("Invalid input: x must have at least one row and one column.")
  }

  if (M <= 0 || N <= 0) {
    stop("Invalid input: M and N must be positive integers.")
  }

  rep_rows <- rep(seq_len(rows), times = M)
  rep_cols <- rep(seq_len(cols), each = N)

  return(x[rep_rows, rep_cols, drop = FALSE])
}

#' Get the classes of variables in a tibble
#'
#' This function takes a tibble \code{x} and returns a character vector containing the classes of each variable in the tibble.
#'
#' @param x A tibble to get the variable classes from
#'
#' @return A character vector containing the classes of each variable in the tibble
#'
#' @export
classes <- function(x) {
  if (!is.data.frame(x) && !is.tibble(x)) {
    stop("Invalid input: x must be a tibble or dataframe.")
  }

  return(sapply(x, class))
}

#' Scale the numeric variables in a tibble
#'
#' This function takes a tibble \code{x} and scales the numeric variables. By default, the variables are centered and scaled to have unit variance. The function returns a tibble with the same attributes as the input tibble.
#'
#' @param x A tibble to scale
#' @param center Logical. Should the variables be centered? (default = TRUE)
#' @param scale Logical. Should the variables be scaled to have unit variance? (default = TRUE)
#'
#' @return A tibble with the numeric variables centered and/or scaled
#'
#' @export
df_scale <- function(x, center = TRUE, scale = TRUE) {
  if (!inherits(x, "tbl_df")) {
    stop("Input must be a tibble")
  }

  num_vars <- sapply(x, is.numeric)

  if (!any(num_vars)) {
    return(x)
  }

  scaled_vars <- scale(x[, num_vars], center = center, scale = scale)

  new_x <- x

  new_x[, num_vars] <- scaled_vars

  return(new_x)
}

#' Calculate the log-likelihood of a sample under the normal distribution
#'
#' This function calculates the log-likelihood of a sample \code{x} under the normal distribution, with mean \code{mean} and standard deviation \code{sd}.
#'
#' @param x A numeric vector of data
#' @param mean The mean of the normal distribution
#' @param sd The standard deviation of the normal distribution
#'
#' @return The log-likelihood of \code{x} under the normal distribution
#'
#' @export
log_likelihood_norm <- function(x, mean, sd) {
  log_likelihood <- sum(dnorm(x, mean, sd, log = TRUE))
  return(log_likelihood)
}

#' Calculate the log-likelihood of a sample under the uniform distribution
#'
#' This function calculates the log-likelihood of a sample \code{x} under the uniform distribution, with minimum value \code{min} and maximum value \code{max}.
#'
#' @param x A numeric vector of data
#' @param min The minimum value of the uniform distribution
#' @param max The maximum value of the uniform distribution
#'
#' @return The log-likelihood of \code{x} under the uniform distribution
#'
#' @export
log_likelihood_unif <- function(x, min, max) {
  log_likelihood <- sum(dunif(x, min, max, log = TRUE))
  return(log_likelihood)
}

#' Calculate the log-likelihood of a sample under the chi-squared distribution
#'
#' This function calculates the log-likelihood of a sample \code{x} under the chi-squared distribution, with degrees of freedom \code{df}.
#'
#' @param x A numeric vector of data
#' @param df The degrees of freedom of the chi-squared distribution
#'
#' @return The log-likelihood of \code{x} under the chi-squared distribution
#'
#' @export
log_likelihood_chisq <- function(x, df) {
  log_likelihood <- sum(dchisq(x, df, log = TRUE))
  return(log_likelihood)
}

#' Calculate the log-likelihood of a sample under the F distribution
#'
#' This function calculates the log-likelihood of a sample \code{x} under the F distribution, with degrees of freedom \code{df1} and \code{df2}.
#'
#' @param x A numeric vector of data
#' @param df1 The degrees of freedom of the numerator of the F distribution
#' @param df2 The degrees of freedom of the denominator of the F distribution
#'
#' @return The log-likelihood of \code{x} under the F distribution
#'
#' @export
log_likelihood_f <- function(x, df1, df2) {
  log_likelihood <- sum(df(x, df1, df2, log = TRUE))
  return(log_likelihood)
}

#' Calculate the log-likelihood of a sample under the t distribution
#'
#' This function calculates the log-likelihood of a sample \code{x} under the t distribution, with degrees of freedom
#' @param x A numeric vector of data
#' @param df The degrees of freedom of the t distribution
#'
#' @return The log-likelihood of \code{x} under the t distribution
#'
#' @export
log_likelihood_t <- function(x, df) {
  log_likelihood <- sum(dt(x, df, log = TRUE))
  return(log_likelihood)
}

#' Calculate sensitivity of a binary classification model
#'
#' This function calculates the sensitivity of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The sensitivity of the binary classification model
#'
#' @export
sensitivity <- function(pred, truth) {
  TP <- sum(pred == 1 & truth == 1)
  FN <- sum(pred == 0 & truth == 1)
  sensitivity <- TP / (TP + FN)
  return(sensitivity)
}

#' Calculate specificity of a binary classification model
#'
#' This function calculates the specificity of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The specificity of the binary classification model
#'
#' @export
specificity <- function(pred, truth) {
  TN <- sum(pred == 0 & truth == 0)
  FP <- sum(pred == 1 & truth == 0)
  specificity <- TN / (TN + FP)
  return(specificity)
}

#' Calculate precision of a binary classification model
#'
#' This function calculates the precision of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The precision of the binary classification model
#'
#' @export
precision <- function(pred, truth) {
  TP <- sum(pred == 1 & truth == 1)
  FP <- sum(pred == 1 & truth == 0)
  precision <- TP / (TP + FP)
  return(precision)
}

#' Calculate recall of a binary classification model
#'
#' This function calculates the recall (sensitivity) of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The recall (sensitivity) of the binary classification model
#'
#' @export
recall <- function(pred, truth) {
  return(sensitivity(pred, truth)) # Recall is the same as sensitivity
}

#' Calculate accuracy of a binary classification model
#'
#' This function calculates the accuracy of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The accuracy of the binary classification model
#'
#' @export
accuracy <- function(pred, truth) {
  correct <- sum(pred == truth)
  accuracy <- correct / length(truth)
  return(accuracy)
}

#' Calculate F1 score of a binary classification model
#'
#' This function calculates the F1 score of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The F1 score of the binary classification model
#'
#' @export
f1 <- function(pred, truth) {
  TP <- sum(pred == 1 & truth == 1)
  FP <- sum(pred == 1 & truth == 0)
  FN <- sum(pred == 0 & truth == 1)
  precision_val <- TP / (TP + FP)
  recall_val <- TP / (TP + FN)
  f1 <- 2 * (precision_val * recall_val) / (precision_val + recall_val)
  return(f1)
}

#' Calculate the minimum n per group needed for a two-sample t-test
#'
#' This function calculates the minimum sample size per group needed for a two-sample t-test, given the expected effect size \code{d} and the desired statistical power \code{power}.
#'
#' @param d The expected Cohen's d
#' @param power The desired statistical power (default 0.8)
#'
#' @return The minimum sample size per group needed for a two-sample t-test
#'
#' @export
minimum_n_per_group <- function(d, power = 0.8) {
  alpha <- 0.05
  n <- pwr.t.test(d = d, power = power, sig.level = alpha, type = "two.sample")$n
  return(ceiling(n))
}

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Calculate R-squared between predicted and true values
#'
#' This function calculates the R-squared statistic between a vector of predicted values \code{pred} and a vector of true values \code{truth}.
#'
#' @param pred A vector of predicted values
#' @param truth A vector of true values
#'
#' @return The R-squared statistic between \code{pred} and \code{truth}
#'
#' @export
r2 <- function(pred, truth) {
  SS_residual <- sum((truth - pred)^2)
  SS_total <- sum((truth - mean(truth))^2)
  r_squared <- 1 - SS_residual / SS_total
  return(r_squared)
}

#' Calculate adjusted R-squared between predicted and true values
#'
#' This function calculates the adjusted R-squared statistic between a vector of predicted values \code{pred} and a vector of true values \code{truth}, given the number of model parameters \code{n_p}.
#'
#' @param pred A vector of predicted values
#' @param truth A vector of true values
#' @param n_p The number of model parameters, excluding the intercept
#'
#' @return The adjusted R-squared statistic between \code{pred} and \code{truth}
#'
#' @export
adj_R2 <- function(pred, truth, n_p) {
  n <- length(truth)
  SS_residual <- sum((truth - pred)^2)
  SS_total <- sum((truth - mean(truth))^2)
  r_squared <- 1 - SS_residual / SS_total
  adj_r_squared <- 1 - ((1 - r_squared) * (n -  1)) / (n - n_p - 1)
  return(adj_r_squared)
}
