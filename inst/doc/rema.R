## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(Rdpack)

## ---- include = FALSE---------------------------------------------------------
trt.events <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 5, 0, 2, 2, 26)
trt.total <- c(121, 59, 358, 454, 238, 62, 186, 87, 369, 257, 798, 978, 605, 
                277, 157, 650)
ctrl.events <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 17, 35, 2, 5, 1, 12)
ctrl.total <- c(118, 58, 164, 216, 268, 59, 97, 94, 386, 109, 804, 996, 608, 
                198, 50, 220)
mid.p <- TRUE
alpha <- 0.05
arguments <- list(as.name("rema"), "rema.obj" = as.name("antibiotics.rema"))
TE <- 0.5708234
CI <-c(0.1947373, 1.5984156)
CMLE <- "CMLE"
pval <- 0.2555818
test.stat <- c(27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 
               43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53)
norm.probs <- c(6.716051e-12, 6.243005e-10, 4.605874e-08, 1.178779e-06, 
                1.810190e-05, 1.900967e-04, 1.361220e-03, 6.418245e-03, 
                2.290489e-02, 5.697199e-02, 1.160254e-01, 1.694432e-01, 
                2.084688e-01, 1.724891e-01, 1.360043e-01, 6.433921e-02, 
                3.335175e-02, 8.452742e-03, 3.156814e-03, 3.155256e-04, 
                8.395551e-05, 2.683688e-06, 7.058532e-07, 6.078317e-09, 
                1.895755e-09, 3.216023e-12, 5.180658e-14)
dist <- data.frame(test.stat, norm.probs)
tstat <- 37

antibiotics.rema <- list("trt.events" = trt.events, 
                         "trt.total" = trt.total, 
                         "ctrl.events" = ctrl.events, 
                         "ctrl.total" = ctrl.total,
                         "mid.p" = mid.p,
                         "alpha" = alpha,
                         "arguments"= arguments, 
                         "TE" = TE,
                         "CI" = CI,
                         "method"= CMLE,
                         "pval" = pval,
                         "dist" = dist,
                         "tstat" = tstat)

class(antibiotics.rema) <- c("rema", "list")

## ----load package, eval=TRUE--------------------------------------------------
library(rema)

## ----Data introduction/background, eval=FALSE---------------------------------
#  rema(trt.events = NULL, trt.total = NULL, ctrl.events = NULL, ctrl.total = NULL,
#       rema.obj, mid.p = TRUE, distr = TRUE, one.sided.p = FALSE, alpha = 0.05)

## ----Data introduction/vectors------------------------------------------------
anti.events <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 5, 0, 2, 2, 26)

anti.total <- c(121, 59, 358, 454, 238, 62, 186, 87, 369, 257, 798, 978, 605, 
                277, 157, 650)

plac.events <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 17, 35, 2, 5, 1, 12)

plac.total <- c(118, 58, 164, 216, 268, 59, 97, 94, 386, 109, 804, 996, 608, 
                198, 50, 220)

## ----Example with default values, eval=FALSE----------------------------------
#  antibiotics.rema <- rema(anti.events, anti.total, plac.events, plac.total)
#  summary(antibiotics.rema)
#  # Call:
#  # rema(trt.events = anti.events, trt.total = anti.total, ctrl.events = plac.events,
#  #     ctrl.total = plac.total)
#  #
#  #         OR           95%-CI p-value
#  #     0.5708 [0.1947; 1.5984]  0.2556
#  #
#  # Details on meta-analytical method:
#  # - Rare event, heterogeneous meta-analysis method
#  # - Two-sided p-value returned (mid.p = TRUE)
#  # - Conditional Maximum Likelihood Estimate (CMLE) used when computing the odds ratio

## ----Example with distr, eval=TRUE--------------------------------------------
antibiotics.rema$dist
antibiotics.rema$tstat

## ----fig.height = 5, fig.width = 5, fig.align = "center"----------------------
plot(antibiotics.rema)

## ----Example with one.sided.p, eval=TRUE--------------------------------------
antibiotics.rema.one.pval <- rema(antibiotics.rema, one.sided.p = TRUE)
antibiotics.rema.one.pval$pval.one.sided

## ----Example with mid.p, eval=TRUE--------------------------------------------
rema(antibiotics.rema, mid.p = FALSE)

## ----Example with alpha, eval=TRUE--------------------------------------------
rema(antibiotics.rema, alpha = 0.1)

