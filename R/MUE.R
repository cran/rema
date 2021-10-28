# Median unbiased estimate by bracketing the root using normalized probabilities

MUE <- function(prob, u, t.obs) {

  f.prob <- prob[u == t.obs] # observed probability
  f.u <- t.obs # observed test statistic

  # if only 1 value of the test statistic in the distribution,
  # then estimation doesn't make a lot of sense
  if (length(u) == 1) {
    return(NA)
  } else {
    # otherwise, use MUE via bracketing the root
    return(FINDROOT(f.prob, f.u, prob, u, rhs = 0.5))
  }
}
