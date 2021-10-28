test_that("throws when vectors are not double or integer", {
  te <- c(NA, NA)
  tt <- c(TRUE, TRUE)
  ce <- c(NA, NA)
  ct <- c(FALSE, FALSE)
  expect_error(rema(te, tt, ce, ct), vector.type.err)

  te <- c("Hi", "Hi")
  tt <- c("Hey", "Hey")
  ce <- c("Hello", "Hello")
  ct <- c("Yo", "Yo")
  expect_error(rema(te, tt, ce, ct), vector.type.err)

  te <- c(1, "Hey")
  tt <- c(2, "Hi")
  ce <- c(3, "Hello")
  ct <- c(4, "Yo")
  expect_error(rema(te, tt, ce, ct), vector.type.err)
})
