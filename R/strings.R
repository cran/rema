# Error strings
min.len.err <- "Input vectors must be of length greater than one."
diff.len.err <- "All input vectors must be of the same length."
vector.type.err <- "All input vectors must be of type double or integer."
vector.val.err <- "Input vectors may not contain Inf, NaN, or NA values."
whole.num.err <- "All input vectors must contain whole numbers. No continuity correction should be applied."
neg.num.err <- "All input vectors must contain non-negative values."
total.size.err <- "Each element in *.events may not exceed the respective element in *.total."
logical.err <- "mid.p, distr, one.sided.pval arguments must be TRUE/FALSE."
alpha.err <- "alpha must satisfy: 0 < alpha < 1."
rema.obj.err <- "rema.obj must be of class 'rema'."
rema.obj.dist.err <- "rema.obj must contain the permutational distribution of the test statistic."
rema.obj.dist.one.val.err <- "There is only one value in the permutational distribution of the test statistic (the observed value), so this method is unable to produce a combined odds ratio and confidence interval."
rema.obj.vectors.err <- "rema.obj must contain the four necessary input vectors with correct labels (trt.events, trt.total, ctrl.events, ctrl.total)"

# Warning strings
obs.freq.warn <- "WARNING: This method is for rare events and may run slowly if the events are more common (roughly, greater than 5% frequency AND large sample sizes)"
test.stat.warn <- "WARNING: Only one value in the permutational distribution (the observed test statistic)."
