test_that("throws when a rema object is given without all input vectors", {
  rema.obj.bad.vectors <- list()
  rema.obj.bad.vectors$dist <- c(1, 2)
  class(rema.obj.bad.vectors) <- "rema"
  expect_error(rema(rema.obj = rema.obj.bad.vectors, rema.obj.vector.err))
})
