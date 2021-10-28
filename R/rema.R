#' A permutation-based approach for heterogeneous meta-analyses of rare events
#'
#' \code{rema} (rare event meta-analysis) performs a permutation-based
#' meta-analysis for heterogeneous, rare event data.
#'
#' Conventional meta-analysis approaches tend to perform poorly for
#' heterogeneous, rare event data. \code{rema} implements a permutation-based
#' approach for binary meta-analyses of 2x2 tables, founded on conditional
#' logistic regression, that provides more reliable statistical tests when
#' heterogeneity is observed in rare event data 
#' \insertCite{Zabriskie2021}{rema}. To adjust for the effect of heterogeneity, 
#' this method conditions on the sufficient statistic of a proxy for the 
#' heterogeneity effect as opposed to estimating the heterogeneity variance. 
#' While this results in the model not strictly falling under the random-effects 
#' framework, it is akin to a random-effects approach in that it assumes 
#' differences in variability due to treatment. Further, this method does not 
#' rely on large-sample approximations or continuity corrections for rare event 
#' data.
#'
#' This method uses the permutational distribution of the test statistic instead
#' of asymptotic approximations for inference. The number of observed events
#' drives the computation complexity for creating this permutational
#' distribution. Accordingly, for this method to be computationally feasible, it
#' should only be applied to meta-analyses with a relatively low number of
#' observed events. To create this permutational distribution, a network
#' algorithm, based on the work of \insertCite{Mehta1992;textual}{rema} and
#' \insertCite{Corcoran2001;textual}{rema}, is employed using C++ and integrated
#' into the package.
#'
#' \insertNoCite{Normand1999}{rema}
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \href{https://www.rdocumentation.org/packages/meta/versions/4.15-1/topics/metabin}{metabin} and
#' \href{https://www.rdocumentation.org/packages/metafor/versions/2.4-0/topics/metafor-package}{metafor}
#' for more traditional meta-analysis methods for combining odds ratios.
#'
#' @param trt.events Numeric vector containing the number of observed events in
#'    the treatment group from each independent study.
#' @param trt.total Numeric vector containing the total number of observations
#'    (events plus non-events) in the treatment group from each independent study.
#' @param ctrl.events Numeric vector containing the number of observed events in
#'    the control group from each independent study.
#' @param ctrl.total Numeric vector containing the total number of observations
#'    (events plus non-events) in the control group from each independent study.
#' @param rema.obj An object of class \code{rema} that contains the permutational distribution
#'    of the test statistic. If this argument is supplied, the trt.events, trt.total,
#'    ctrl.events, and ctrl.total vectors are not required.
#' @param mid.p A logical indicating if the p-values and confidence interval are
#'    adjusted using the mid-p correction to reduce conservatism. If \code{TRUE}
#'    (default), the mid-p p-values and mid-p confidence interval will be
#'    reported.
#' @param distr A logical indicating if the permutational distribution of the
#'    test statistic is reported. If \code{TRUE} (default), the distribution
#'    is returned.
#' @param one.sided.p A logical indicating if the one-sided p-value (in the
#'    direction of the treatment effect) is computed and returned. If
#'    \code{FALSE} (default), the one-sided p-value is not provided.
#' @param alpha A number between 0 and 1 to create a (1 - alpha)% confidence
#'    interval for the treatment effect.
#'
#' @return An object of \code{\link{class}} "\code{rema}" with corresponding
#'   \code{print}, \code{summary}, and \code{plot} (if \code{distr == TRUE})
#'   functions. The object is a list containing the following elements:
#'   \tabular{ll}{
#'     \code{ trt.events} \tab As defined above. \cr
#'     \code{trt.total} \tab As defined above. \cr
#'     \code{ctrl.events} \tab As defined above. \cr
#'     \code{ctrl.total} \tab As defined above. \cr
#'     \code{mid.p} \tab As defined above. \cr
#'     \code{alpha} \tab As defined above. \cr
#'     \code{arguments} \tab A string of the arguments passed into \code{rema}. \cr
#'     \code{TE} \tab The estimated overall treatment effect (odds ratio). \cr
#'     \code{CI} \tab A vector containing the estimated lower and upper
#'                            bounds of a (1 - alpha)% \cr 
#'                            \tab confidence interval of the
#'                            overall treatment effect (odds ratio). \cr
#'     \code{method} \tab A string specifying the method used to compute
#'                          the odds ratio and its associated \cr
#'                          \tab confidence 
#'                          interval. Either the conditional maximum likelihood 
#'                          estimate (CMLE) \cr
#'                          \tab or the median unbiased estimate (MUE)
#'                          will be used. \cr
#'     \code{pval} \tab The two-sided p-value for the overall treatment effect. \cr
#'     \code{pval.one.sided} \tab The one-sided p-value (in the direction of
#'                                  the treatment effect) for the
#'                                   overall \cr 
#'                                   \tab treatment effect (if
#'                                   \code{one.sided.p == TRUE}). \cr
#'     \code{dist} \tab A data frame containing the permutational
#'                                distribution of the test statistic \cr
#'                                \tab (if
#'                                \code{distr == TRUE}). \cr
#'     \code{tstat} \tab The observed value of the test statistic (if 
#'                         \code{distr == TRUE}). \cr
#'    }
#'
#' @export
#'
#' @examples
#' # Lidocaine data set from Normand (1999)
#' rema(trt.events = c(2, 4, 6, 7, 7, 11),
#'      trt.total = c(39, 44, 107, 103, 110, 154),
#'      ctrl.events = c(1, 4, 4, 5, 3, 4),
#'      ctrl.total = c(43, 44, 110, 100, 106, 146),
#'      mid.p = FALSE,
#'      distr = FALSE,
#'      one.sided.p = TRUE)
#'
#' # Example using a rema object as the input once the permutational
#' # distribution is obtained
#' my.rema.object <- rema(trt.events = c(1, 2, 0, 2, 0, 0, 1, 0, 1, 0),
#'                        trt.total = c(30, 14, 30, 49, 38, 11, 31, 13, 49, 23),
#'                        ctrl.events = c(4, 3, 3, 0, 4, 4, 3, 2, 3, 4),
#'                        ctrl.total = c(15, 26, 42, 24, 40, 26, 47, 24, 27, 26),
#'                        mid.p = FALSE,
#'                        distr = TRUE,
#'                        one.sided.p = FALSE)
#' rema(rema.obj = my.rema.object,
#'      mid.p = TRUE)
#'
#' \dontrun{
#' # Vectors of non-whole numbers (such as after applying a continuity correction)
#' rema(trt.events = c(0.5, 0.5, 1, 3),
#'      trt.total = c(2, 5, 4, 12),
#'      ctrl.events = c(2.5, 4.5, 6, 7),
#'      ctrl.total = c(7, 9, 11, 12))
#'
#' # Vectors with greater observed events than total observations
#' rema(trt.events = c(11, 13, 7, 10),
#'      trt.total = c(10, 12, 5, 7),
#'      ctrl.events = c(22, 25, 32, 26),
#'      ctrl.total = c(20, 20, 30, 25))
#' }
#'

rema <- function(trt.events = NULL, trt.total = NULL, ctrl.events = NULL,
                 ctrl.total = NULL, rema.obj, mid.p = TRUE, distr = TRUE,
                 one.sided.p = FALSE, alpha = 0.05) {
  # Check arguments
  # if only trt.events vector is given check to see if it was rema.obj
  if ((!is.null(trt.events))
     && is.null(trt.total)
     && is.null(ctrl.events)
     && is.null(ctrl.total)
     && ("rema" %in% class(trt.events))) {
    rema.obj <- trt.events
  }
  use.rema.obj <- check.inputs(trt.events, trt.total, ctrl.events, ctrl.total,
                               rema.obj, mid.p, distr, one.sided.p, alpha)



  # pass in to c++ code or use given distribution
  if (use.rema.obj) {
    distribution <- rema.obj$dist
    num.studies <- length(rema.obj$trt.events)
    treated <- c(rep(0, num.studies), rep(1, num.studies))
    ctrl.trt.total <- c(rema.obj$ctrl.total, rema.obj$trt.total)
    ctrl.trt.events <- c(rema.obj$ctrl.events, rema.obj$trt.events)
    trt.events = rema.obj$trt.events
    trt.total = rema.obj$trt.total
    ctrl.events = rema.obj$ctrl.events
    ctrl.total = rema.obj$ctrl.total
  } else {
    # prep data
    num.studies <- length(trt.events)
    samp.size <- sum(ctrl.total, trt.total)
    treated <- c(rep(0, num.studies), rep(1, num.studies))
    events <- sum(ctrl.events, trt.events)

    ob.stat <- sum(trt.events)
    obs.corr.0 <- sum(ctrl.events * (ctrl.total - ctrl.events))
    obs.corr.1 <- sum(trt.events * (trt.total - trt.events))

    p.val <- 0

    ctrl.trt.total <- c(ctrl.total, trt.total)
    ctrl.trt.events <- c(ctrl.events, trt.events)

    distribution <- trstatWrapper(num.studies * 2, ctrl.trt.total, treated,
                                  samp.size, events, num.studies, ob.stat,
                                  obs.corr.0, obs.corr.1, p.val)
    # check distribution
    if (nrow(distribution) <= 1) {
      stop(rema.obj.dist.one.val.err)
    }

    # remove unnorm.probs column
    distribution$unnorm.probs <- NULL
  }

  # pass distribution into inference function
  input.data <- data.frame(trt = treated,
                           event = ctrl.trt.events,
                           total = ctrl.trt.total)

  final.results <- rema.inference(distribution, input.data, mid.p, one.sided.p,
                                  alpha)

  final.results$arguments <- match.call()

  return.results <- list("trt.events" = trt.events,
                         "trt.total" = trt.total,
                         "ctrl.events" = ctrl.events,
                         "ctrl.total" = ctrl.total,
                         "mid.p" = mid.p)

  return.results$alpha <- final.results$alpha
  return.results$arguments <- final.results$arguments

  if (all(!is.na(final.results))) {
    # Exponentiate to return to original scale
    return.results$TE <- exp(final.results$TE)
    return.results$CI <- exp(final.results$CI)
  }

  return.results$method <- final.results$method
  return.results$pval <- final.results$pval
  return.results$pval.one.sided <- final.results$pval.one.sided

  # Include distribution and tstat if distr == TRUE
  if (distr) {
    return.results$dist <- distribution
    return.results$tstat <- final.results$tstat
  }

  class(return.results) <- c("rema", "list")
  check.outputs(return.results)
  return(return.results)
}


check.outputs <- function(remaObj) {
  if (!is.null(remaObj$distribution)) {
    distr <- remaObj$distribution
    if (nrow(distr) == 1) {
      message(test.stat.warn)
    }
  }
}


check.inputs <- function(trt.events, trt.total, ctrl.events, ctrl.total,
                         rema.obj, mid.p, distr, one.sided.p, alpha) {
  use.rema.obj <- FALSE

  if (!missing(rema.obj)) {
    # check to see if rema object is actually of class rema
    if (!("rema" %in% class(rema.obj))) {
      stop(rema.obj.err)
    }
    # check to see if rema object has distribution
    if (is.null(rema.obj$dist)) {
      stop(rema.obj.dist.err)
    }
    if (is.null(rema.obj$trt.events)
        || is.null(rema.obj$trt.total)
        || is.null(rema.obj$ctrl.events)
        || is.null(rema.obj$ctrl.total)) {
      stop(rema.obj.vectors.err)
    }
    trt.events <- rema.obj$trt.events
    trt.total <- rema.obj$trt.total
    ctrl.events <- rema.obj$ctrl.events
    ctrl.total <- rema.obj$ctrl.total
    use.rema.obj <- TRUE
  }
  check.vector.inputs(trt.events, trt.total, ctrl.events, ctrl.total)
  if (!is.logical(mid.p)
      || is.na(mid.p)) {
    stop(logical.err)
  }
  if (!is.logical(distr)
      || is.na(distr)) {
    stop(logical.err)
  }
  if (!is.logical(one.sided.p)
      || is.na(one.sided.p)) {
    stop(logical.err)
  }
  if (!is.double(alpha)
      || is.na(alpha)
      || 0 >= alpha
      || 1 <= alpha) {
    stop(alpha.err)
  }
  return(use.rema.obj)
}

check.vector.inputs <- function(trt.events, trt.total, ctrl.events, ctrl.total) {
  if (length(trt.events) <= 1) {
    stop(min.len.err)
  }
  if (length(trt.events) != length(trt.total)
      || length(trt.total) != length(ctrl.events)
      || length(ctrl.events) != length(ctrl.total)) {
    stop(diff.len.err)
  }
  if ((!is.double(trt.events) && !is.integer(trt.events))
      || (!is.double(trt.total) && !is.integer(trt.total))
      || (!is.double(ctrl.events) && !is.integer(ctrl.events))
      || (!is.double(ctrl.total) && !is.integer(ctrl.total))) {
    stop(vector.type.err)
  }
  if (any(is.infinite(trt.events))
      || any(is.nan(trt.events))
      || any(is.na(trt.events))
      || any(is.infinite(trt.total))
      || any(is.nan(trt.total))
      || any(is.na(trt.total))
      || any(is.infinite(ctrl.events))
      || any(is.nan(ctrl.events))
      || any(is.na(ctrl.events))
      || any(is.infinite(ctrl.total))
      || any(is.nan(ctrl.total))
      || any(is.na(ctrl.total))) {
    stop(vector.val.err)
  }
  if (any(trt.events %% 1 != 0)
      || any(trt.total %% 1 != 0)
      || any(ctrl.events %% 1 != 0)
      || any(ctrl.total %% 1 != 0)) {
    stop(whole.num.err)
  }
  if (any(trt.events < 0)
      || any(trt.total < 0)
      || any(ctrl.events < 0)
      || any(ctrl.total < 0)) {
    stop(neg.num.err)
  }
  if (any(trt.events > trt.total)
      || any(ctrl.events > ctrl.total)) {
    stop(total.size.err)
  }
  if (((mean(trt.events / trt.total) > 0.05
          || mean(ctrl.events / ctrl.total) > 0.05)
        & (any(trt.total > 100)
          || any(ctrl.total > 100)))
      || (mean(trt.events / trt.total) > 0.15
        || mean(ctrl.events / ctrl.total) > 0.15)) {
    message(obs.freq.warn)
  }
}
