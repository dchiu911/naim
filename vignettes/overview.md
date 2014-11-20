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
