test_that("throws when number of events exceed total observations", {
  te <- c(2, 2)
  tt <- c(1, 1)
  ce <- c(3, 3)
  ct <- c(0, 0)
  # This is a subset of the error message total.size.err to comply with the
  # regex expectations of the string search
  expect_error(rema(te, tt, ce, ct),
               " may not exceed the respective element in ")
})
