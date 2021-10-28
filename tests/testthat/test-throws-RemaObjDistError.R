test_that("throws when a rema object is given without distribution", {
  rema.obj.no.dist <- c(1, 2)
  class(rema.obj.no.dist) <- "rema"
  expect_error(rema(rema.obj = rema.obj.no.dist, rema.obj.dist.err))
})
