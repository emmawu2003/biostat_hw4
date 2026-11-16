
# matrix_lm

<!-- badges: start -->
<!-- badges: end -->

The goal of matrix_lm is to fits ordinary least squares regression models using QR
decomposition. It provides `matrix_lm()`, which takes a model formula
and a data frame. The function omits rows with missing values in the
model variables and requires a design matrix with full column rank.

## Installation

Install the development version from GitHub:

``` r
# If needed, first run install.packages("remotes").
remotes::install_github("emmawu2003/biostat_hw4")
```

## Example

Use the built-in `mtcars` dataset to fit a linear model for fuel economy
(`mpg`) as a function of vehicle weight (`wt`).

```r
library(simplelm)

fit <- matrix_lm(mpg ~ wt, data = mtcars)
fit$coefficients
```

The output contains an intercept and a slope for `wt`. The slope is the
estimated change in miles per gallon for each additional 1,000 pounds
of vehicle weight.

The fitted model returns a list with three components:

- `coefficients` — the estimated regression coefficients.
- `Q` — the matrix with orthonormal columns from the QR decomposition.
- `R` — the upper triangular matrix from the QR decomposition.

For more details on the arguments and return values, run `?matrix_lm`
after loading the package.

## Comparison with R's `lm()`

Both `matrix_lm()` and `stats::lm()` use QR decomposition to fit
ordinary least squares models. `matrix_lm()` constructs the QR
decomposition using modified Gram–Schmidt and returns the estimated
coefficients and the `Q` and `R` matrices directly. `stats::lm()` returns
an `"lm"` object with fitted values and residuals and supports methods
such as `summary()` and `predict()`.

`matrix_lm()` requires a full-rank design matrix. `stats::lm()` allows
rank-deficient fits by default.

For the `mtcars` example above, the two fits should give the same
coefficients up to numerical rounding:

```r
base_fit <- stats::lm(mpg ~ wt, data = mtcars)
all.equal(fit$coefficients, coef(base_fit))
```

See the [R documentation for `lm()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/lm.html)
for more details.
