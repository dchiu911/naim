#' Newton-Raphson method for computing MLE under logistic regression model
#'
#' For a simple logistic regression model with an intercept and one slope parameter, this function returns their respective MLE's. Iterations stop once the improvement in successive steps falls under a certain tolerance level.
#'
#' @details It is known that in logsitic regression, we cannot obtain closed form solutions for the parameters. This function uses the Newton-Raphson method for the (slope, paramter) vector. Our initial value is (0, 0), an intuitive choice that corresponds to equally likely outcomes. We can view \code{y} as the number of successes out of \code{n} trials for each point.
#' testing
#'
#' @param x vector of covariates (predictor)
#' @param y vector of responses, each one distributed as a Binomial random variable
#' @param n number of trials for each sample
#' @param tol tolerance level governing when to stop the iterations. Here we use the larger absolute error for either parameter and compare it to \code{tol}. Iterations stop when the error is less than \code{tol}.
#' @param verbose Logical; if \code{TRUE} then function prints useful feedback to the console, and if \code{FALSE} there is no printing. Defaults to \code{TRUE}
#' @return A vector of length two, the first being the MLE of the intercept and the second being the MLE of the slope.
#' @keywords misc
#'
NR_logit <- function(x, y, n, tol, verbose = T){
	X <- cbind(1, x)
	d <- c(0, 0)
	if(verbose = T){
		print(paste("Initial Value:", "alpha =", d[1], "beta =", d[2]))
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
		if(verbose = T){
			print(paste("Iteration", it, ":", "alpha =", round(d_new[1], 6),
									", beta =", round(d_new[2], 6)))
		}
	}
	alpha_hat <- round(d_new[1], 6)
	beta_hat <- round(d_new[2], 6)
	if(verbose = T){
		print(paste("MLE are alpha =", alpha_hat, "and beta =", beta_hat))
	}
	return(d_new)
}
