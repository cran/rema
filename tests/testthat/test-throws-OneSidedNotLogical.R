test_that("throws when one.sided.p is not TRUE/FALSE", {
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), one.sided.p = 1), logical.err)
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), one.sided.p = 1.1), logical.err)
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), one.sided.p = "Hi"), logical.err)
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), one.sided.p = NA), logical.err)
})
