test_that("throws when vectors are smaller than minimum length", {
  te <- c()
  tt <- c()
  ce <- c()
  ct <- c()
  expect_error(rema(te, tt, ce, ct), min.len.err)

  te <- c(1)
  tt <- c(2)
  ce <- c(3)
  ct <- c(4)
  expect_error(rema(te, tt, ce, ct), min.len.err)
})
