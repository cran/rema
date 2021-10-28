rema.inference <- function(distribution, input.data, mid.p, one.sided.p,
                           alpha) {
  # observed test statistic from the data
  t.obs <- as.numeric(input.data$trt %*% input.data$event)

  #########################################################################
  #                      compute two-sided p-value
  #########################################################################

  obs.prob <- distribution$norm.probs[distribution$test.stat == t.obs]
  if (mid.p) {
    pval <- sum(distribution$norm.probs[distribution$norm.probs < obs.prob]) +
      (0.5 * distribution$norm.probs[distribution$norm.probs == obs.prob])
  } else {
    pval <- sum(distribution$norm.probs[distribution$norm.probs <= obs.prob])
  }

  #########################################################################
  #                      compute one-sided p-values
  #########################################################################
  # split distribution based on observed test statistic
  if (one.sided.p) {

    # report the one-sided p-value in the direction of the observed test stat
    if (length(distribution[distribution$test.stat < t.obs, 1]) <
        length(distribution[distribution$test.stat > t.obs, 1])) {
      part.dist <- distribution[distribution$test.stat <= t.obs, ]
    } else {
      part.dist <- distribution[distribution$test.stat >= t.obs, ]
    }

    if (mid.p) {
      one.sided.p.value <-
        sum(part.dist$norm.probs[part.dist$norm.probs < obs.prob]) +
        (0.5 * part.dist$norm.probs[part.dist$norm.probs == obs.prob])
    } else {
      one.sided.p.value <-
        sum(part.dist$norm.probs[part.dist$norm.probs <= obs.prob])
    }
  }

  #########################################################################
  #                        save parts of distribution
  #########################################################################
  u <- distribution$test.stat
  prob <- ifelse(distribution$norm.probs == 0,
                 1.0e-16,
                 distribution$norm.probs)

  #########################################################################
  #        note if t.obs is at either extreme of the distribution
  #########################################################################
  t.obs.at.extremes <- ifelse(t.obs == max(u) || t.obs == min(u),
                              TRUE, FALSE)
  t.obs.approx.at.extremes <- ifelse(sum(prob[u > t.obs]) < 1.0e-15 ||
                                       sum(prob[u < t.obs]) < 1.0e-15,
                                     TRUE, FALSE)

  #########################################################################
  #      if t.obs is at either extreme of the distribution, use MUE
  #########################################################################
  if (t.obs.at.extremes || t.obs.approx.at.extremes) {

    beta.hats <- MUE(prob = prob,
                     u = u,
                     t.obs = t.obs)
    ci <- CI(prob = prob,
             u = u,
             t.obs = t.obs,
             beta = beta.hats,
             alpha = alpha,
             midp = mid.p)
    conf.low <- ci[1]
    conf.up <- ci[2]

    method <- "MUE"

    #########################################################################
    #        otherwise (if t.obs is NOT at an extreme), use CMLE
    #########################################################################
  } else {

    # scale down test statistics to start at 0
    test.obs <- t.obs - min(distribution$test.stat)
    u <- distribution$test.stat - min(distribution$test.stat)

    beta.hats <- CMLE(prob = prob,
                      u = u,
                      t.obs = test.obs,
                      t.or.alpha = test.obs,
                      ci = FALSE)
    conf.low <- CMLE(prob = prob,
                     u = u,
                     t.obs = test.obs,
                     t.or.alpha = alpha/2,
                     ci = TRUE,
                     dir = "lower",
                     midp = mid.p)
    conf.up <- CMLE(prob = prob,
                    u = u,
                    t.obs = test.obs,
                    t.or.alpha = alpha/2,
                    ci = TRUE,
                    dir = "upper",
                    midp = mid.p)

    method <- "CMLE"
  }


  # all on log scale
  final.results <- list(TE = beta.hats,
                        CI = c(conf.low, conf.up),
                        tstat = t.obs,
                        pval = pval,
                        alpha = alpha,
                        method = method)

  if (one.sided.p) {
    final.results$pval.one.sided <- one.sided.p.value
  }

  return(final.results)
}
