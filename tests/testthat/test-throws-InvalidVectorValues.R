test_that("throws when input vectors contain invalid values", {
  expect_error(rema(c(1, 3), c(2, Inf), c(0, 1), c(1, 2), mid.p = 1), vector.val.err)
  expect_error(rema(c(1, 3), c(2, NaN), c(0, 1), c(1, 2), mid.p = 1), vector.val.err)
  expect_error(rema(c(1, 3), c(2, NA), c(0, 1), c(1, 2), mid.p = 1), vector.val.err)
})
