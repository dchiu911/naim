## Tests for `NR_logit`
# Correct input vectors
x <- runif(15, 6, 10)
n <- sample(1:15, replace = TRUE)
y <- rbinom(15, n, 0.6)

# Check data.frame return
test_that("Always returns data.frame", {
	expect_is(NR_logit(x, y, n), "data.frame")
})

# Check error for different cases of invalid input vectors
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
