test_that("DAPSm with fixed weight for toyData2", {
  data('toyData2')
  toyData2$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
                              data = toyData2)$fitted.values
  daps <- DAPSest(toyData2, out.col = 2, trt.col = 1, caliper = 0.3,
                  weight = 0.7, coords.columns = c(4, 5),
                  pairsRet = TRUE, cov.cols = 6:9, cutoff = 0.1,
                  coord_dist = TRUE, caliper_type = 'DAPS',
                  matching_algorithm = 'greedy')
  expect_equal(daps$weight, 0.7)
})


test_that("DAPSm with fast optimal for toyData2", {
  data('toyData2')
  toyData2$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
                              data = toyData2)$fitted.values
  daps <- DAPSest(toyData2, out.col = 2, trt.col = 1, caliper = 0.3,
                  weight = 'optimal', coords.columns = c(4, 5), quiet = TRUE,
                  pairsRet = TRUE, cov.cols = 6:9, cutoff = 0.15,
                  coord_dist = TRUE, caliper_type = 'DAPS', w_tol = 0.05,
                  matching_algorithm = 'greedy')
  expect_equal(daps$weight, 0.515625)
  expect_equal(as.numeric(abs(daps$est - 1.2868) < 0.001), 1)
})


test_that("DAPSm with extensive optimal for toyData2", {
  data('toyData2')
  toyData2$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
                              data = toyData2)$fitted.values
  bal <- CalcDAPSWeightBalance(toyData2, weights = seq(0, 1, length.out = 40),
                               cov.cols = 6:9, trt.col = 1,
                               coords.columns = c(4, 5), caliper = 0.3,
                               matching_algorithm = 'greedy')
  daps <- DAPSchoiceModel(toyData, trt.col = 1, balance = bal$balance,
                          cutoff = 0.15, pairs = bal$pairs,
                          weights = seq(0, 1, length.out = 40))
  expect_equal(abs(daps$weight - 0.2307692) < 0.001, TRUE)
  expect_equal(daps$num_match, 55)
  print(abs(daps$est - 2.53958) < 0.1)
  expect_equal(as.numeric(abs(daps$est - 2.5) < 1), 1)
})
