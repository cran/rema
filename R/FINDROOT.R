# Bracketing the root using normalized probabilities

FINDROOT <- function(f.prob, f.u, prob, u, rhs, lower = -1, upper = 1) {

  prob <- log(prob) # normalized probabilities on log scale
  f.prob <- log(f.prob) # normalized probabilities on log scale

  exact.cond.prob.fun <- function(f.prob, f.u, prob, u, beta) {
    snew <- prob[1] + u[1] * beta
    if (length(u) > 1) {
      for (i in 2:length(u)) {
        # add log probabilities
        mxt <- max(snew, prob[i] + u[i] * beta)
        snew <- mxt + log(exp(snew - mxt)
                          + exp(prob[i] + u[i] * beta - mxt))
      }
    }

    f.snew <- f.prob[1] + f.u[1] * beta
    if (length(f.u) > 1) {
      for (i in 2:length(f.u)) {
        f.mxt <- max(f.snew, f.prob[i] + f.u[i] * beta)
        f.snew <- f.mxt + log(exp(f.snew - f.mxt) +
                              exp(f.prob[i] + f.u[i] * beta - f.mxt))
      }
    }
    return(exp(f.snew - snew))
  }


  a <- lower # Starting value for lower bound
  b <- upper # Starting value for upper bound
  f.of.a <- exact.cond.prob.fun(f.prob, f.u, prob, u, a) - rhs
  f.of.b <- exact.cond.prob.fun(f.prob, f.u, prob, u, b) - rhs

  # change b (increase magnitude and flip sign)
  # until a and b bracket a root
  # i.e. until f.of.a and f.of.b are of opposite signs
  while (sign(f.of.a) == sign(f.of.b)) {
    b <- abs(b) + 1
    f.of.b <- exact.cond.prob.fun(f.prob, f.u, prob, u, b) - rhs
    if (sign(f.of.a) == sign(f.of.b)) {
      b <- -b
      f.of.b <- exact.cond.prob.fun(f.prob, f.u, prob, u, b) - rhs
    }
  }

  # now that we have a and b that bracket a root,
  # narrow down to the root's value
  stop <- 0
  num.steps <- 0

  while (stop == 0 && num.steps < 10000) {
    root <- (f.of.b * a - f.of.a * b) / (f.of.b - f.of.a)
    f.of.ck <- exact.cond.prob.fun(f.prob, f.u, prob, u, root) - rhs
    num.steps <- num.steps + 1
    if (f.of.ck == 0) {
      stop <- 1
    } else {
      if (sign(f.of.a) == sign(f.of.ck)) {
        a <- root
        f.of.a <- f.of.ck
      } else {
        b <- root
        f.of.b <- f.of.ck
      }
      if ((abs(a - b) < 1e-10)) {
        stop <- 1
      }
    }
  }

  return(root)
}
