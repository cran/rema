test_that("throws when a non rema object is given", {
  expect_error(rema(rema.obj = 1, rema.obj.err))
  expect_error(rema(rema.obj = NA, rema.obj.err))
  expect_error(rema(rema.obj = FALSE, rema.obj.err))
  expect_error(rema(rema.obj = NaN, rema.obj.err))
  expect_error(rema(rema.obj = "rema", rema.obj.err))
})
