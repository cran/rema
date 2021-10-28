test_that("throws when distr is not TRUE/FALSE", {
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), distr = 1), logical.err)
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), distr = 1.1), logical.err)
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), distr = "Hi"), logical.err)
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), distr = NA), logical.err)
})
