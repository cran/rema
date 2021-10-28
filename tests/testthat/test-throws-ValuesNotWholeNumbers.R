test_that("throws when vector values aren't whole numbers", {
  te <- c(1.1, 1.1)
  tt <- c(2.5, 2.5)
  ce <- c(3.2, 3.5)
  ct <- c(6.8, 6.8)
  expect_error(rema(te, tt, ce, ct), whole.num.err)
})
