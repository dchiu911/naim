#' Newton-Raphson method for computing MLE under logistic regression model
#'
#' For a simple logistic regression model with an intercept and one slope parameter, this function returns their respective MLE's. Iterations stop once the improvement in successive steps falls under a certain tolerance level.
#'
#'
#' We can view \code{y} as the number of successes out of \code{n} trials, with probability of success \code{p} unknown. The linear predictor is related to \code{p} by the log odds ratio: \eqn{log(p/(1-p)) = \alpha + \beta x}.
#' It is known that in logsitic regression, we cannot obtain closed form solutions for \eqn{\alpha} and \eqn{\beta}. This function uses the Newton-Raphson method to estimate these parameters. Our initial value is (\eqn{\alpha_0}, \eqn{\beta_0)} = (0, 0), an intuitive choice that corresponds to the case of equally likely outcomes (p = 0.5).
#'
#' In the function definition, the matrix \code{L} is the system of nonlinear score functions we need to solve, and \code{L_prime} is the derivative of \code{L}. Hence,each step of the iteration computes (\eqn{\alpha_{i+1}}, \eqn{\beta_{i+1}}) = (\eqn{\alpha_{i}}, \eqn{\beta_{i}}) + \code{L_prime}^-1\ * \code{L}.
#'
#' @param x vector of covariates (predictor)
#' @param y vector of responses, each one distributed as a Binomial random variable
#' @param n vector of number of trials for each observation.
#' @param tol tolerance level governing when to stop the iterations. Here we use the larger absolute error for either parameter and compare it to \code{tol}. Iterations stop when the error is less than \code{tol}. Defaults to \code{1e-6}.
#' @param verbose Logical; if \code{TRUE} then function prints iterative feedback to the console, and if \code{FALSE} there is no printing. Defaults to \code{FALSE}.
#' @return A data.frame with two entries in one row, the first being the MLE of the intercept and the second being the MLE of the slope.
#' @export
#' @examples
#' # Generate random covariates, responses, and trial totals
#' set.seed(547)
#' x <- runif(15, 6, 10)
#' n <- sample(1:15, replace = TRUE)
#' y <- rbinom(15, n, 0.6)
#' NR_logit(x, y, n)
#' NR_logit(x, y, n, verbose = TRUE)
#' NR_logit(x, y, n, tol = 1e-10, verbose = TRUE)
NR_logit <- function(x, y, n, tol = 1e-6, verbose = FALSE){
	assertthat::assert_that(check_vectors(x, y, n))

	X <- cbind(1, x)
	d <- c(0, 0)

	if(verbose == TRUE){
		print(paste("Initial Value:", "intercept =", d[1], "slope =", d[2]))
	}
	it = 0
	tolerance = 1

	while(tolerance > tol) {
		p <- exp(X %*% d) / (1 + exp(X %*% d))
		V <- diag(c(n*p*(1 - p)))
		L <- t(X) %*% (y - n * p)
		L_prime <- t(X) %*% V %*% X
		d_new <- d + solve(L_prime) %*% L

		tolerance <- max(abs(d_new[1] - d[1]), abs(d_new[2] - d[2]))
		d <- d_new
		it <- it + 1

		if(verbose == TRUE){
			print(paste("Iteration", it, ":", "intercept =", round(d_new[1], 6),
									", slope =", round(d_new[2], 6)))
		}
	}
	d_new <- data.frame(intercept = d_new[1], slope = d_new[2])

	intercept_hat <- round(d_new[1], 6)
	slope_hat <- round(d_new[2], 6)

	if(verbose == TRUE){
		print(paste("MLE are intercept =", intercept_hat, "and slope =", slope_hat))
	}
	return(d_new)
}
