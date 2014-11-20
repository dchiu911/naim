#' EM algorithm for estimating frequency of A, B, and O blood alleles
#'
#' Given the phenotypic frequencies for the four blood types (A, B, AB, O), uses the EM algorithm to obtain maximum likelihood estimates for the frequency of the three blood alleles (A, B, O) in a population.
#'
#' The genotypic frequencies are usually not known for an entire population, because people with blood type A can have a genotype of A/A or A/O, and people with blood type B can have a genotype of B/B or B/O. For blood types AB and O, there is a one-to-one mapping to the genotypes A/B and O/O however. Hence, in the EM setting, the six possible genotypes are viewed as latent variables. Using traditional MLE, closed form solutions cannot obtained.
#'
#' We assume that the population is in Hardy-Weinberg Equilibrium. That is, for genotype of form X/X, P(X/X) = P(X)^2 and for genotype of form X/Y, P(X/Y) = 2 * P(X) * P(Y).
#'
#' An initial value of (1/3, 1/3, 1/3) is reasonable as it corresponds to the situation of equal allele frequency. In the function definition, the iterations are updated in the calls to \code{pA_new}, \code{pB_new}, \code{pO_new}. The equations can be verified by going through the Expectation-Maximization algorithm.
#'
#' @param A number of people with blood type A
#' @param B number of people with blood type B
#' @param AB number of people with blood type AB
#' @param O number of people with blood type O
#' @param tol tolerance level governing when to stop the iterations. Here we use the large absolute error for any parameter and compare it to \code{tol}. Iterations stop when the error is less than \code{tol}. Defaults to \code{1e-6}.
#' @param verbose Logical; if \code{TRUE} then function prints iterative feedback to the console, and if \code{FALSE} there is no printing. Defaults to \code{FALSE}.
#' @return A data.frame with three entries in one row, giving the estimated allele frequencies of A, B, and O respectively.
#' @export
#' @examples
#' # Population of 100, with equal phenotypic frequencies
#' A <- 25; B <- 25; AB <- 25; O <- 25
#' EM_blood(A, B, AB, O)
#'
#' # Population of 1000, with different phenotypic frequencies
#' A <- 500; B <- 200; AB <- 50; O <- 250
#' EM_blood(A, B, AB, O)
EM_blood <- function(A, B, AB, O, tol = 1e-6, verbose = FALSE){
	assertthat::assert_that(all(!is.na(c(A, B, AB, O))))

	n <- A + B + AB + O
	p <- round(rep(1/3, 3), 6)
	pA <- p[1]; pB <- p[2]; pO <- p[3]

	if(verbose == TRUE){
		print(paste("Initial Probs:", "pA =", pA, "pB =", pB, "pO = ", pO))
	}
	it = 0
	tolerance = 1

	while(tolerance > tol){
		A.A <- A * pA^2/(pA^2 + 2*pA*pO)
		B.B <- B * pB^2/(pB^2 + 2*pB*pO)
		O.A <- A * 2*pA*pO/(pA^2 + 2*pA*pO)
		O.B <- B * 2*pB*pO/(pB^2 + 2*pB*pO)

		pA_new <- (2*A.A + AB + O.A)/(2*n)
		pB_new <- (AB + 2*B.B + O.B)/(2*n)
		pO_new <- (2*O + O.A + O.B)/(2*n)

		tolerance = max(abs(pA_new - pA), abs(pB_new - pB), abs(pO_new - pO))
		pA <- pA_new
		pB <- pB_new
		pO <- pO_new
		it <- it + 1
		if(verbose == TRUE){
			print(paste("Iteration", it, ":", "pA =", round(pA_new, 6),
									"pB =", round(pB_new, 6), "pO =", round(pO_new, 6)))
		}
	}
	pA_hat <- round(pA_new, 6)
	pB_hat <- round(pB_new, 6)
	pO_hat <- round(pO_new, 6)
	if(verbose == TRUE){
		print(paste("MLE are pA =", pA_hat, "pB =", pB_hat, "pO =", pO_hat))
	}
	return(data.frame(pA_hat = pA_hat, pB_hat = pB_hat, pO_hat = pO_hat))
}
