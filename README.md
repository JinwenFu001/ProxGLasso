# ProxGLasso

`ProxGLasso` is an R package for sparse inverse covariance estimation in Gaussian graphical models. It focuses on the graphical lasso problem and provides proximal optimization methods for estimating a sparse precision matrix.

The package includes a backtracking proximal gradient method, a self-concordant proximal gradient method, and proximal Newton methods based on both primal and dual updates. It is intended for settings where standard first-order methods may be difficult to apply directly because of the log-determinant term in the objective.

## Installation

You can install the development version from GitHub with:

```r
install.packages("remotes")
remotes::install_github("JinwenFu001/ProxGLasso", build_vignettes = TRUE)
```

Then load the package with:

~~~r
library(ProxGLasso)
~~~

## Problem overview

In a Gaussian graphical model, the goal is to estimate a precision matrix whose off-diagonal nonzero entries encode conditional dependence relationships among variables. The graphical lasso does this by combining a Gaussian log-likelihood term with an elementwise `l1` penalty on the off-diagonal entries of the precision matrix.

`ProxGLasso` implements proximal methods for this optimization problem. In particular, it includes methods designed for the graphical lasso setting where the smooth part of the objective involves a log-determinant term and is not globally Lipschitz differentiable.

## Main functions

`Generate_AR1_pair()` generates Gaussian samples under an AR1 covariance or precision structure.

`BTProx_glasso()` fits the graphical lasso using a proximal gradient method with backtracking line search.

`SCProx_glasso()` fits the graphical lasso using a self-concordant proximal gradient method.

`Newton_Prox_primal_glasso()` and `Newton_Prox_dual_glasso()` fit the graphical lasso using proximal Newton methods based on primal and dual formulations.

## Minimal example

The example below generates Gaussian data, computes the sample covariance matrix, and fits the graphical lasso using several methods implemented in the package.

```r
library(ProxGLasso)

set.seed(42)

n <- 10
p <- 5

samples <- Generate_AR1_pair(n, p, rho = 0.7, cov_AR1 = TRUE)
S <- cov(samples$X) / n * (n - 1)

BT_est <- BTProx_glasso(S, n, lambda = 0.1)
SC_est <- SCProx_glasso(S, n, lambda = 0.1)

SC_Newton_pri_est <- Newton_Prox_primal_glasso(S, n, lambda = 0.1)
SC_Newton_dual_est <- Newton_Prox_dual_glasso(S, n, lambda = 0.1)

BT_est
SC_est
SC_Newton_pri_est
SC_Newton_dual_est
```

## Documentation

For a longer introduction and additional methodological details, see the package vignette:

```r
browseVignettes("ProxGLasso")
```

## License

This package is distributed under the license specified in the `LICENSE.md` file.