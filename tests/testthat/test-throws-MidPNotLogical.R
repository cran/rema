test_that("throws when mid.p is not logical", {
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), mid.p = 1), logical.err)
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), mid.p = 1.1), logical.err)
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), mid.p = "Hi"), logical.err)
  expect_error(rema(c(1, 1), c(2, 2), c(0, 0), c(1, 1), mid.p = NA), logical.err)
})
