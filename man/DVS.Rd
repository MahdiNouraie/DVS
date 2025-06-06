% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DVS.R
\docType{package}
\name{DVS}
\alias{-package}
\alias{DVS}
\alias{getStability}
\title{DVS: Decorrelation for Variable Selection}
\usage{
DVS(x, y, B)
}
\arguments{
\item{x}{A numeric matrix of predictors.}

\item{y}{A numeric vector of response values.}

\item{B}{An integer specifying the number of sub-samples in the stability selection.}
}
\value{
`lambda.stable` and its associated stability value, along with variables whose selection frequencies exceed 0.5 and their corresponding frequencies. If `lambda.stable` is not attainable, the function uses `lambda.stable.1sd` instead.
}
\description{
This package contains a main function: `DVS`.
The `DVS` function first decorrelates variables and then applies stability selection with Lasso to find relevant variables.
The function prints `lambda.stable` and its associated stability value, along with variables whose selection frequencies exceed 0.5 and their corresponding frequencies. If `lambda.stable` is not attainable, the function uses `lambda.stable.1sd` instead.
}
\examples{
set.seed(123)
n <- 100 # Number of observations
rho <- 0.8 # Correlation coefficient for the predictors
x1 <- matrix(rnorm(n * 3), ncol = 3) # First 3 independent predictors
x2 <- rho * x1[, rep(1:3, length.out = 7)] + sqrt(1 - rho^2) * matrix(rnorm(n * 7), ncol = 7) # Make next 7 predictors correlated with x1
x <- cbind(x1, x2) # Combine independent and correlated predictors
colnames(x) <- paste0("X", 1:10) # Assign column names
beta <- c(1, 2, 3, rep(0, 7)) # Create regression coefficients vector
y <- x \%*\% beta + rnorm(n) # Generate response variable with some noise
B <- 10 # Number of sub-samples for stability selection
DVS(x, y, B)  # Example usage of the DVS function
#output
$lambda.stable
[1] 0.2798887
$stability
[1] 0.7781636
Variable Selection_Frequency
1       X3                   1
2       X2                   1
3       X1                   1

}
\references{
Joudah, I., Muller, S., & Zhu, H. (2025). Air-HOLP: adaptive regularized feature screening for high dimensional correlated data. Statistics and Computing, 35(3), 63.

Nouraie, M., & Muller, S. (2024). On the Selection Stability of Stability Selection and Its Applications. arXiv preprint arXiv:2411.09097.

Nogueira, S., Sechidis, K., & Brown, G. (2018). On the stability of feature selection algorithms. Journal of Machine Learning Research, 18(174), 1-54.

Meinshausen, N., & Bühlmann, P. (2010). Stability selection. Journal of the Royal Statistical Society Series B: Statistical Methodology, 72(4), 417-473.

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society Series B: Statistical Methodology, 58(1), 267-288.
}
\author{
Mahdi Nouraie (mahdinouraie20@gmail.com)
}
