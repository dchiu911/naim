naim
=======

[![Build Status](https://travis-ci.org/dchiu911/naim.svg?branch=master)](https://travis-ci.org/dchiu911/naim)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dchiu911/naim?branch=master&svg=true)](https://ci.appveyor.com/project/dchiu911/naim)
[![codecov](https://codecov.io/gh/dchiu911/naim/branch/master/graph/badge.svg)](https://codecov.io/gh/dchiu911/naim)

A package that uses <b>n</b>umerical <b>a</b>nalysis <b>i</b>terative <b>m</b>ethods to find maximum likelihood estimates when closed form solutions cannot be obtained.

### Installation
Install `naim` from GitHub:


```r
devtools::install_github("dchiu911/naim")
```
_Note that you need the `devtools` package to do this._

### Usage
First, load the package the usual way:


```r
library(naim)
```

At the current version, `naim` has two functions: `NR_logit()` and `ER_blood()`.

Given vectors `x`, `n`, and `y`, for the covariates, trials, and successes in a logistic regression setting, the MLE for the intercept and slope computed using the Newton-Raphson method are:


```r
set.seed(1)
x <- rnorm(100, mean = 3, sd = 0.2)
n <- sample(1:100, replace = TRUE)
y <- rbinom(100, size = n, prob = 0.6)
NR_logit(x, y, n)
```

```
##   intercept      slope
## 1  0.639122 -0.0694504
```

Given the number of people with blood type A, B, AB, and O, the frequency of the blood alleles A, B, and O are:


```r
A <- 10
B <- 20
AB <- 30
O <- 40
EM_blood(A, B, AB, O)
```

```
##    pA_hat   pB_hat   pO_hat
## 1 0.20833 0.270615 0.521055
```

Both functions return `data.frames` with the pertinent information in labelled columns.

### Vignette
Further details can be found in the vignette for this package. To view it, run this command in the R console:


```r
browseVignettes("naim")
```

Once there, you can click on the `HTML` link to get a nice introduction to `naim`. The help files can be accessed with `?NR_logit` and `?EM_blood` for details on function arguments and some examples.

Alternatively, there is a markdown version available in this repository, found [here](https://github.com/dchiu911/naim/blob/master/vignettes/overview.md).
