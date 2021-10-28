#' Plots the distribution of the test statistic results from a \code{rema} object
#'
#' Plots the distribution of the test statistic results from a \code{rema} object 
#' with the probabilities on the y-axis and the test statistics on the x-axis.
#'
#' The observed value of the test statistic is marked with an asterisk. The dark 
#' gray bars mark the probabilities that contribute to the two-sided p-value. 
#' The light gray bars highlight probabilities that are greater than the 
#' probability of the observed test statistic and consequently do not contribute 
#' to the two-sided p-value. If \code{dist} is set to \code{FALSE} in the 
#' \code{rema} object, the plot function will return an error.
#'
#' @param x a \code{rema} object holding distribution results from the 
#'    \code{rema} function
#' @param ... other arguments that can be specified in plot
#'
#' @return A barplot of the distribution of the test statistic.
#' 
#' @export
#'
plot.rema <- function(x, ...) {
  if (!(is.null(x$dist) || is.null(x$tstat))) {
    colors <- ifelse(x$dist$norm.probs <= x$dist$norm.probs[x$dist$test.stat == x$tstat],
                     "gray47", "gray82")
    bar.text <- ifelse(x$dist$test.stat == x$tstat, "*", "")
    remaPlot <- graphics::barplot(x$dist$norm.probs,
                                  names.arg = x$dist$test.stat,
                                  col = colors,
                                  ylim = c(0, ifelse(max(x$dist$norm.probs) < 0.1,
                                                     max(x$dist$norm.probs) + 0.01,
                                                     max(x$dist$norm.probs) + 0.05)),
                                  ylab = "Probability", xlab = "Test Statistic")
    graphics::text(remaPlot, x$dist$norm.probs + 0.01, bar.text, font = 2)
  } else {
    stop("Plot is not available. Please rerun the rema function with distr == TRUE, and then try again.")
  }
}


#' @export
print.rema <- function(x, ...) {
  # Convert to data frame for printing purposes
  x.df <- data.frame(beta.hat = x$TE,
                     conf.up = x$CI[2],
                     conf.low = x$CI[1],
                     pval = x$pval)

  if (all(!is.na(x.df))) {
    # Round to 4 decimal places
    x.df <- (round(x.df, digits = 4))
  }

  # Format (1 - alpha)% confidence interval and append to data frame
  x.df$CI <- paste(c("[", x.df[[2]]), c(x.df[[3]], "]"),
                   sep = "", collapse = "; ")

  # Edit row names and column names
  rownames(x.df) <- "   "
  ci.column.name <- paste(c((1 - x$alpha) * 100, "%-CI"), collapse = "")
  colnames(x.df) <- c("OR", "conf.up", "conf.low", "p-value", ci.column.name)

  # Print information about the data
  cat("Call:\n")
  print(x$arguments)
  cat("\n")

  # Print OR, (1 - alpha)% CI, and p-value
  print.data.frame(x.df[, c(1, 5, 4)])
  cat("\n")

  # Print details about the method
  cat("Details on meta-analytical method:\n")
  cat("- Rare event, heterogeneous meta-analysis method\n")
  if (x$mid.p) {
    cat("- Two-sided p-value returned (mid.p = TRUE)\n")
  } else {
    cat("- Two-sided p-value returned (mid.p = FALSE)\n")
  }

  if (x$method == "MUE") {
    cat("- Median Unbiased Estimate (MUE) used when computing the odds ratio\n")
  } else {
    cat("- Conditional Maximum Likelihood Estimate (CMLE) used when computing the odds ratio\n")
  }
}


#' @export
summary.rema <- function(object, ...) {
  # Convert to data frame for printing purposes
  object.df <- data.frame(beta.hat = object$TE,
                          conf.up = object$CI[2],
                          conf.low = object$CI[1],
                          pval = object$pval)

  if (all(!is.na(object.df))) {
    # Round to 4 decimal places
    object.df <- (round(object.df, digits = 4))
  }

  # Format (1 - alpha)% confidence interval and append to data frame
  object.df$CI <- paste(c("[", object.df[[2]]), c(object.df[[3]], "]"),
                        sep = "", collapse = "; ")

  # Edit row names and column names
  rownames(object.df) <- "   "
  ci.column.name <- paste(c((1 - object$alpha) * 100, "%-CI"), collapse = "")
  colnames(object.df) <- c("OR", "conf.up", "conf.low", "p-value", ci.column.name)

  # Print information about the data
  cat("Call:\n")
  print(object$arguments)
  cat("\n")

  # Print OR, (1 - alpha)% CI, and p-value
  print.data.frame(object.df[, c(1, 5, 4)])
  cat("\n")

  # Print details about the method
  cat("Details on meta-analytical method:\n")
  cat("- Rare event, heterogeneous meta-analysis method\n")
  if (object$mid.p) {
    cat("- Two-sided p-value returned (mid.p = TRUE)\n")
  } else {
    cat("- Two-sided p-value returned (mid.p = FALSE)\n")
  }

  if (object$method == "MUE") {
    cat("- Median Unbiased Estimate (MUE) used when computing the odds ratio\n")
  } else {
    cat("- Conditional Maximum Likelihood Estimate (CMLE) used when computing the odds ratio\n")
  }
}