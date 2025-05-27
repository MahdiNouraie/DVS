# DVS: Decorrelation for Variable Selection

**DVS** is an R package designed for stable variable selection in the presence of correlated predictors using Lasso within the stability selection framework.

The methodology is based on the paper:  
**"Stability Selection via Variable Decorrelation" (2025) — Nouraie et al.**

## Key Features

- Performs variable selection using Lasso under stability selection.
- Applies a variable decorrelation step prior to selection to enhance stability.
- Returns selected variables along with their selection frequencies and stable regularisation parameter.

## Installation


DVS depends on the R packages `glmnet` and `cmna`, which will be automatically installed when you install `DVS`.

You can install and load the `DVS` package using the following commands in R:

```r
# Install 'devtools' if not already installed
if (!require("devtools")) {
  install.packages("devtools")
}

# Install the DVS package from GitHub
devtools::install_github("MahdiNouraie/DVS")

# Load the package
library(DVS)
```

## Example Usage

```r
set.seed(123)

# Simulate data with correlated predictors
n <- 100                         # Number of observations
rho <- 0.8                       # Desired correlation between predictors
x1 <- matrix(rnorm(n * 3), ncol = 3)  # First 3 independent predictors
x2 <- rho * x1[, rep(1:3, length.out = 7)] + 
      sqrt(1 - rho^2) * matrix(rnorm(n * 7), ncol = 7)  # Correlated predictors
x <- cbind(x1, x2)               # Combine all predictors
colnames(x) <- paste0("X", 1:10)

# Generate response
beta <- c(1, 2, 3, rep(0, 7))    # True coefficients
y <- x %*% beta + rnorm(n)       # Response with noise

# Apply DVS
B <- 10                          # Number of subsamples
DVS(x, y, B)
```

### Example Output

```
$lambda.stable
[1] 0.2798887

$stability
[1] 0.7781636

$Variable_Selection_Frequency
  Variable Frequency
1       X3         1
2       X2         1
3       X1         1
```

## Acknowledgements

`DVS` includes adapted code from the following sources, which are appropriately cited in the code with comments:
- [JMLR2018 Supplementary Code](https://github.com/nogueirs/JMLR2018)
- [Air-HOLP Repository](https://github.com/Logic314/Air-HOLP)
- [StackOverflow – Gram-Schmidt in R](https://stackoverflow.com/questions/15584221/gram-schmidt-with-r)

## License

This package is released under the MIT License.
