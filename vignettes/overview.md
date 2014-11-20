### Motivation

The usual method of finding maximum likelihood estimates involves
deriving the log likelihood function with respect to each parameter we
want to estimate, then setting this equation to 0. But there are time
when this is not possible. For example, in the logistic regression
scenario, closed form solutions cannot be obtained. Thus we need to
resort to some iterative method in order to compute these estimates. In
this way, we usually say the MLE is found when the improvement after
each step is smaller than some tolerance level.

This package currently implements the Newton-Raphson method for simple
logistic regression, and the EM algorithm for estimating the frequency
of blood alleles.

### Usage

First load the package via:

    library(naim)

The `NR_logit` function find estimated parameters under simple logistic
regression, where closed form solutions are not obtainable. Using a
method like Newton-Raphson, we supply an initial feasible estimate, and
iteratively improve the estimate.

All we need is to supply it with **3 required arguments**: covariates,
number of trials, and successes. Details found in the documentation via
`?NR_logit`. Let's see an example:

    # Covariates
    x <- rnorm(100, mean = 3, sd = 0.2)

    # Number of trials
    n <- sample(1:100, replace = TRUE)

    # Successes
    y <- rbinom(100, size = n, prob = 0.6)

    # Data set
    head(data.frame(x, n, y))

    ##          x  n  y
    ## 1 3.221225 82 52
    ## 2 3.095797 29 15
    ## 3 2.880289 87 50
    ## 4 2.961991 91 56
    ## 5 2.874805 68 38
    ## 6 2.623090 84 53

The default call to `NR_logit` returns a data.frame with the intercept
and slope MLE.

    NR_logit(x, y, n)

    ##   intercept      slope
    ## 1 0.3565655 0.01299847

Suppose we wanted to see the iterative process, i.e. the incremental
improvements in the estimated parameters. This can be done by setting
the argument `verbose = TRUE`.

    NR_logit(x, y, n, verbose = TRUE)

    ## [1] "Initial Value: intercept = 0 slope = 0"
    ## [1] "Iteration 1 : intercept = 0.352972 , slope = 0.012503"
    ## [1] "Iteration 2 : intercept = 0.356565 , slope = 0.012998"
    ## [1] "Iteration 3 : intercept = 0.356565 , slope = 0.012998"
    ## [1] "MLE are intercept = 0.356565 and slope = 0.012998"

    ##   intercept      slope
    ## 1 0.3565655 0.01299847

It shows the initial value used (which is always (0, 0)), the updated
estimates at each iteration, and the final MLE.

The stopping criterion occurs when the largest absolute error of either
parameter falls below a certain tolerance level. By default the
tolerance is set to **1e-6**. If we wanted to increase the accuracy of
our estimate, we can manually set a lower tolerance by adding the
argument `tolerance = 1e-10`, for example.

    NR_logit(x, y, n, tol = 1e-15, verbose = TRUE)

    ## [1] "Initial Value: intercept = 0 slope = 0"
    ## [1] "Iteration 1 : intercept = 0.352972 , slope = 0.012503"
    ## [1] "Iteration 2 : intercept = 0.356565 , slope = 0.012998"
    ## [1] "Iteration 3 : intercept = 0.356565 , slope = 0.012998"
    ## [1] "Iteration 4 : intercept = 0.356565 , slope = 0.012998"
    ## [1] "Iteration 5 : intercept = 0.356565 , slope = 0.012998"
    ## [1] "Iteration 6 : intercept = 0.356565 , slope = 0.012998"
    ## [1] "Iteration 7 : intercept = 0.356565 , slope = 0.012998"
    ## [1] "MLE are intercept = 0.356565 and slope = 0.012998"

    ##   intercept      slope
    ## 1 0.3565655 0.01299847

Notice that the Newton-Raphson takes 12 iterations to converge now
because the tolerance has been set to an extremely small value. In the
context of this example, there is no significant difference from our
previous estimate.

The second function of this package uses the EM algorithm to estimate
the frequency of blood alleles (A, B, or O) in a population given their
phenotypic frequencies. An iterative method is needed because closed
form solutions cannot be obtained. The allele frequencies are
intertwined within the six possible blood genotypes, forming the latent
variables in our problem.

Using `EM_blood`, you can quickly obtain the three MLE by providing the
total number of people who have A, B, AB, and O. The arguments `tol` and
`verbose` are also implemented and function the exact same way as in
`NR_logit`.

Let's work through a simple example. Suppose the phenotypic frequencies
are as given below:

    A <- 80; B <- 45; AB <- 13; O <- 100

Then the allele frequencies are:

    EM_blood(A, B, AB, O)

    ##     pA_hat   pB_hat   pO_hat
    ## 1 0.219678 0.130473 0.649849

The prevalence of the B antigen is the lowest in this population whereas
having no antigen has the highest probability. *Note that the sum of the
probabilities is always 1 (with some rounding error).*

Again, we can expect the iterative process using `verbose = TRUE`:

    EM_blood(A, B, AB, O, verbose = TRUE)

    ## [1] "Initial Probs: pA = 0.333333 pB = 0.333333 pO =  0.333333"
    ## [1] "Iteration 1 : pA = 0.251401 pB = 0.153361 pO = 0.595238"
    ## [1] "Iteration 2 : pA = 0.224682 pB = 0.132638 pO = 0.642681"
    ## [1] "Iteration 3 : pA = 0.220385 pB = 0.130692 pO = 0.648923"
    ## [1] "Iteration 4 : pA = 0.219775 pB = 0.130498 pO = 0.649728"
    ## [1] "Iteration 5 : pA = 0.219691 pB = 0.130476 pO = 0.649833"
    ## [1] "Iteration 6 : pA = 0.21968 pB = 0.130474 pO = 0.649846"
    ## [1] "Iteration 7 : pA = 0.219678 pB = 0.130473 pO = 0.649848"
    ## [1] "Iteration 8 : pA = 0.219678 pB = 0.130473 pO = 0.649849"
    ## [1] "MLE are pA = 0.219678 pB = 0.130473 pO = 0.649849"

    ##     pA_hat   pB_hat   pO_hat
    ## 1 0.219678 0.130473 0.649849
