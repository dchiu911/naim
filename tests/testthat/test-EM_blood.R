
context("EM_blood")

A <- 80
B <- 150
AB <- 30
O <- 200

test_that("Always returns data.frame", {
	expect_is(EM_blood(A, B, AB, O), "data.frame")
})

test_that("Sum of three estimated frequencies is 1", {
	expect_equal(sum(EM_blood(A, B, AB, O)), 1)
})

test_that("Error if any argument has NA", {
	expect_error(EM_blood(NA, B, AB, O))
	expect_error(EM_blood(A, NA, AB, O))
	expect_error(EM_blood(A, B, NA, O))
	expect_error(EM_blood(A, B, AB, NA))
})

test_that("Feedback prints to console", {
	expect_output(EM_blood(A, B, AB, O, verbose = TRUE))
})
