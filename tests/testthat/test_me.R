## Tests for `NR_logit`
x <- runif(15, 6, 10)
n <- sample(1:15, replace = TRUE)
y <- rbinom(15, n, 0.6)

test_that("Always returns data.frame", {
	expect_is(NR_logit(x, y, n), "data.frame")
})

not_vector <- TRUE
too_short <- x[-1]
has_NA <- c(x, NA)[-1]
y_larger_than_n <- y + n

test_that("Can't operate on invalid input vectors.", {
	expect_error(NR_logit(not_vector, y, n))
	expect_error(NR_logit(too_short, y, n))
	expect_error(NR_logit(has_NA, y, n))
	expect_error(NR_logit(x, y_larger_than_n, n))
})


## Tests for `EM_blood`
A <- 80; B <- 150; AB <- 30; O = 200

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
