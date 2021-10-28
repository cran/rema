# Exact (1 - alpha)% confidence interval based on median unbiased estimation
# with normalized probabilities

CI <- function(prob, u, t.obs, beta, alpha = 0.05, midp) {
  if (is.na(beta)) {
    result <- c(NA, NA)

  } else {
    if (sum(prob[u > t.obs]) < 1.0e-15) {
      prob <- prob[u <= t.obs]
      u <- u[u <= t.obs]
    }
    if (sum(prob[u < t.obs]) < 1.0e-15) {
      prob <- prob[u >= t.obs]
      u <- u[u >= t.obs]
    }

    # Upper tail of test statistic distribution
    # Distribution from t.obs to t.max

    if (midp) {
      f.prob <- c(0.5 * prob[u == t.obs], prob[u > t.obs])
    } else {
      f.prob <- prob[u >= t.obs]
    }

    f.u <- u[u >= t.obs]

    if (t.obs == min(u)) {
      CI.1 <- -Inf
    } else {
      CI.1 <- FINDROOT(f.prob, f.u, prob, u, rhs = 0.5 * alpha,
                       lower = beta - 0.5, upper = beta + 0.5)
    }

    # Lower tail of test statistic distribution
    # Distribution from t.min to t.obs

    if (midp) {
      f.prob <- c(0.5 * prob[u == t.obs], prob[u < t.obs])
    } else {
      f.prob <- prob[u <= t.obs]
    }

    f.u <- u[u <= t.obs]

    if (t.obs == max(u)) {
      CI.2 <- Inf
    } else {
      CI.2 <- FINDROOT(f.prob, f.u, prob, u, rhs = 0.5 * alpha,
                       lower = beta - 0.5, upper = beta + 0.5)
    }

    result <- c(min(CI.1, CI.2), max(CI.1, CI.2))
  }

  return(result)
}
