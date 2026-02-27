library(blockmodels)

expect_model_estimation <- function(model) {
  expect_no_error(model$estimate())
  expect_true(is.numeric(model$ICL))
  expect_gte(length(model$ICL), 1)
  expect_type(which.max(model$ICL), "integer")
}
