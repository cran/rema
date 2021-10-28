# rema

The `rema` (rare event meta-analysis) package implements a permutation-based approach for binary
meta-analyses of 2x2 tables, founded on conditional logistic regression, that 
provides more reliable statistical tests when heterogeneity is observed in rare 
event data (Zabriskie et al. 2021). To adjust for the effect of heterogeneity, this method conditions on
the sufficient statistic of a proxy for the heterogeneity effect as opposed to 
estimating the heterogeneity variance. While this results in the model not 
strictly falling under the random-effects framework, it is akin to a 
random-effects approach in that it assumes differences in variability due to 
treatment. Further, this method does not rely on large-sample approximations or 
continuity corrections for rare event data. 

This method uses the permutational 
distribution of the test statistic instead of asymptotic approximations for inference. The number of observed events
drives the computation complexity for creating this permutational
distribution. Accordingly, for this method to be computationally feasible, it
should only be applied to meta-analyses with a relatively low number of
observed events. To create this permutational distribution, a network
algorithm, based on the work of Mehta et al. (1992) and
Corcoran et al. (2001), is employed using C++ and integrated
into the package.

## Installation

You can install the released version of `rema` from
[CRAN](https://CRAN.R-project.org) with

``` r
install.packages("rema")
```

and the development version from [GitHub](https://github.com/) through
the `devtools` package with

``` r
# install.packages("devtools")
devtools::install_github("ZabStatLab/rema")
```

## Usage

Here is an example using `rema` for a rare event meta-analysis on a small
example data set. For more detailed analyses, please refer to the package vignette.

``` r
library(rema)

trt.events <- c(2, 4, 6, 7, 7, 11)
trt.total <- c(39, 44, 107, 103, 110, 154)
ctrl.events <- c(1, 4, 4, 5, 3, 4)
ctrl.total <- c(43, 44, 110, 100, 106, 146)

rema(trt.events, trt.total, ctrl.events, ctrl.total)                                                                         
#> Call:
#> rema(trt.events = trt.events, trt.total = trt.total, ctrl.events = ctrl.events, 
#>     ctrl.total = ctrl.total)
#> 
#>         OR           95%-CI p-value
#>     0.6457 [0.1512; 3.2015]  0.2423
#> 
#> Details on meta-analytical method:
#> - Rare event, heterogeneous meta-analysis method
#> - Two-sided p-value returned (mid.p = TRUE)
#> - Conditional Maximum Likelihood Estimate (CMLE) used when computing the odds ratio
```

## References

Corcoran C, Ryan L, Senchaudhuri P, Mehta C, Patel N, Molenberghs G (2001). “An Exact Trend
Test for Correlated Binary Data.” _Biometrics_, 57, 941–948, doi:10.1111/j.0006-341x.2001.00941.x.

Mehta CR, Patel N, Senchaudhuri P (1992). “Exact Stratified Linear Rank Tests for Ordered Categorical
and Binary Data.” _Journal of Computational and Graphical Statistics_, 1(1), 21–40, doi:10.2307/1390598.

Zabriskie BN, Corcoran C, Senchaudhuri P (2021). "A Permutation-Based Approach for 
Heterogeneous Meta-Analyses of Rare Events." _Statistics in Medicine_, 40(25), 5587-5604,
doi:10.1002/sim.9142.
