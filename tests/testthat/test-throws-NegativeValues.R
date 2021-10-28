test_that("throws when values are negative", {
  te <- c(-2, 2)
  tt <- c(-1, 1)
  ce <- c(-5, 5)
  ct <- c(-2, 2)
  expect_error(rema(te, tt, ce, ct), neg.num.err)
})
