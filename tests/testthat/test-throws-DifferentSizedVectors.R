test_that("throws when input vectors are not the same size", {
  te <- c(0, 2, 3)
  tt <- c(4, 3, 4)
  ce <- c(3, 4)
  ct <- c(6, 6)
  expect_error(rema(te, tt, ce, ct), diff.len.err)

  tt <- c(1, NULL)
  expect_error(rema(te, tt, ce, ct), diff.len.err)
})
