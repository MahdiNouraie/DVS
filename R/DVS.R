#' @title DVS: Decorrelation for Variable Selection
#' @description
#' This package contains a main function: `DVS`.
#' The `DVS` function first decorrelates predictors using Air-HOLP and Gram-Schmidt orthogonalization, and then applies stability selection with the Lasso to identify relevant variables.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom cmna newton
#' @importFrom stats coef qnorm sd
#' @export
#' @param x A numeric matrix of predictors.
#' @param y A numeric vector of response values.
#' @param B An integer specifying the number of sub-samples in the stability selection.
#' @param Threshold An integer specifying the number of variables retained by
#' Air-HOLP during the initial screening step before decorrelation. The `Threshold` parameter controls the number of variables retained after
#' Air-HOLP screening. It should typically be chosen as a small multiple of the
#' expected number of relevant variables.
#'
#' @return
#' A list containing either:
#' \itemize{
#'   \item \code{lambda.stable} and its associated stability value, if a stable tuning parameter is found;
#'   \item otherwise, \code{lambda.stable.1sd} and its associated stability value.
#' }
#' In both cases, the returned list also contains a data frame of selected variables and their selection frequencies.
#'
#'
#' @examples
#' set.seed(123)
#' n <- 100 # Number of observations
#' rho <- 0.8 # Correlation coefficient for the predictors
#' x1 <- matrix(rnorm(n * 3), ncol = 3) # First 3 independent predictors
#' x2 <- rho * x1[, rep(1:3, length.out = 7)] +
#' sqrt(1 - rho^2) * matrix(rnorm(n * 7), ncol = 7) # Make next 7 predictors correlated with x1
#' x <- cbind(x1, x2) # Combine independent and correlated predictors
#' colnames(x) <- paste0("X", 1:10) # Assign column names
#' beta <- c(1, 2, 3, rep(0, 7)) # Create regression coefficients vector
#' y <- x %*% beta + rnorm(n) # Generate response variable with some noise
#' B <- 10 # Number of sub-samples for stability selection
#' DVS(x, y, B, Threshold = 10)  # Example usage of the DVS function
#'


#' @references
#' Joudah, I., Muller, S., & Zhu, H. (2025). Air-HOLP: adaptive regularized feature screening for high dimensional correlated data. Statistics and Computing, 35(3), 63.
#'
#' Nouraie, M., & Muller, S. (2024). On the Selection Stability of Stability Selection and Its Applications. arXiv preprint arXiv:2411.09097.
#'
#' Nogueira, S., Sechidis, K., & Brown, G. (2018). On the stability of feature selection algorithms. Journal of Machine Learning Research, 18(174), 1-54.
#'
#' Meinshausen, N., & Bühlmann, P. (2010). Stability selection. Journal of the Royal Statistical Society Series B: Statistical Methodology, 72(4), 417-473.
#'
#' Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society Series B: Statistical Methodology, 58(1), 267-288.
#'
#'
#'
#'
DVS <- function(x, y, B, Threshold = 10){
  if (!is.matrix(x) || !is.numeric(x))
    stop("'x' must be a numeric matrix.")

  if (!is.numeric(y) || length(y) != nrow(x))
    stop("'y' must be a numeric vector with length equal to nrow(x).")

  if (!is.numeric(B) || length(B) != 1 || B < 2 || B != as.integer(B))
    stop("'B' must be a single integer >= 2.")

  if (!is.numeric(Threshold) || length(Threshold) != 1 ||
      Threshold < 1 || Threshold > ncol(x) ||
      Threshold != as.integer(Threshold))
    stop("'Threshold' must be a single integer between 1 and ncol(x).")

  x <- scale(x)  # Standardize the predictors
  y <- scale(y, scale = FALSE) # Center the response
  p <- ncol(x) # Number of predictors
  AHOLP <- AirHOLP(x, y, Threshold = Threshold, r0 = 10, adapt = TRUE, iter = 10)
  ranked_features <- AHOLP$index_r  # Ranking of features
  column_order <- order(ranked_features)  # Find the order that sorts ranks from 1 to p

  x_ranked <- x[, column_order]   # Reorder columns of X
  var_names <- colnames(x_ranked) # Store names of predictors
  gram <- gram_schmidtR(x_ranked) # Perform Gram-Schmidt orthogonalization
  x <- gram$Q # Extract orthogonalized matrix
  colnames(x) <- var_names # Set column names for orthogonalized matrix
  cv_lasso <- cv.glmnet(x, y, nfolds = 10, alpha = 1) # Fit Lasso model with 10-fold CV
  candidate_set <- cv_lasso$lambda # Candidate set of lambda values
  S_list <- vector("list", length(candidate_set)) # Initialize a list to store selection matrix for each lambda
  names(S_list) <- paste0("lambda_", seq_along(candidate_set)) # Name the list entries

  for (lambda_idx in seq_along(candidate_set)) { # Stability Selection for each lambda in candidate_set

    lambda <- candidate_set[lambda_idx]  # Current lambda value
    S <- matrix(0, nrow = B, ncol = p)  # Initialize selection matrix for the current lambda
    colnames(S) <- colnames(x) # Set column names of S to predictor names

    for (i in 1:B) {
      # sub-sampling the data (half of the original data without replacement)
      sample_index <- sample(seq_len(nrow(x)), floor(nrow(x) / 2))

      # Prepare the response and predictors
      x_sub <- x[sample_index,] # Predictor variables
      y_sub <-y[sample_index] # Response variable

      # Fit the Lasso model with the current lambda
      lasso_model <- glmnet(x_sub, y_sub, alpha = 1, lambda = lambda)

      # Extract significant predictors (ignoring the intercept, hence [-1])
      significant_predictors <- ifelse(coef(lasso_model) != 0, 1, 0)[-1]

      # Store the significant predictors in matrix S
      S[i, ] <- significant_predictors
    }

    # Store the matrix S for the current lambda in the corresponding list entry
    S_list[[lambda_idx]] <- S
  }
  stability_results <- lapply(S_list, getStability)
  stab_values <- unlist(lapply(stability_results, function(x) x$stability))

  if (all(is.na(stab_values))) {
    stop("All stability estimates are NA.")
  }

  if (max(stab_values, na.rm = TRUE) > 0.75) {
    stable_values <- which(stab_values > 0.75) # Index of stable lambda values
    lambda_stable <- min(candidate_set[stable_values]) # Minimum stable lambda value
    index_of_lambda_stable <- which(candidate_set == lambda_stable) # Index of lambda_stable
    phi_stable <- stab_values[index_of_lambda_stable] # Stability value for lambda_stable
    result <- list(
      lambda.stable = lambda_stable,
      stability = as.numeric(phi_stable)
    )

    Stable_S <- S_list[[index_of_lambda_stable]] # Stable selection matrix for lambda_stable
    col_means <- colMeans(Stable_S)
    selected_cols <- col_means[col_means > 0.5] # Select columns with mean > 0.5
    selected_df <- data.frame(
      Variable = names(selected_cols),
      Selection_Frequency = as.numeric(selected_cols),
      row.names = NULL
    )
    result$selected <- selected_df
    return(result)
  }
  else{
    max_stability <- max(stab_values, na.rm = TRUE) # Find the maximum stability value
    stability_1sd_threshold <- max_stability - sd(stab_values, na.rm = TRUE) # Define the stability threshold as max stability - 1SD
    index_of_stable_1sd <- max(which(stab_values >= stability_1sd_threshold), na.rm = TRUE) # since candidate values are sorted decreasingly, we use max index
    lambda_stable_1sd <- candidate_set[index_of_stable_1sd] # Find the corresponding lambda value
    phi_stable_1sd <- stab_values[index_of_stable_1sd] # Find the corresponding stability value
    result <- list(
      lambda.stable.1sd = lambda_stable_1sd,
      stability = as.numeric(phi_stable_1sd)
    )

    Stable_S_1sd <- S_list[[index_of_stable_1sd]] # Find the corresponding selection matrix
    col_means_1sd <- colMeans(Stable_S_1sd)
    selected_cols_1sd <- col_means_1sd[col_means_1sd > 0.5] # Select columns with mean > 0.5
    selected_df <- data.frame(
      Variable = names(selected_cols_1sd),
      Selection_Frequency = as.numeric(selected_cols_1sd),
      row.names = NULL
    )
    result$selected <- selected_df
    return(result)
  }
}

# Measuring selection stability
#
# Stability measure (2018) from "https://github.com/nogueirs/JMLR2018/blob/master/R/getStability.R"
#
#' @keywords internal
#' @noRd
getStability <- function(X,alpha=0.05) {
  ## the input X is a binary matrix of size M*d where:
  ## M is the number of bootstrap replicates
  ## d is the total number of features
  ## alpha is the level of significance (e.g. if alpha=0.05, we will get 95% confidence intervals)
  ## it's an optional argument and is set to 5% by default
  ### first we compute the stability

  M<-nrow(X)
  d<-ncol(X)
  hatPF<-colMeans(X)
  kbar<-sum(hatPF)
  v_rand=(kbar/d)*(1-kbar/d)
  stability<-1-(M/(M-1))*mean(hatPF*(1-hatPF))/v_rand ## this is the stability estimate

  ## then we compute the variance of the estimate
  ki<-rowSums(X)
  phi_i<-rep(0,M)
  for(i in 1:M){
    phi_i[i]<-(1/v_rand)*((1/d)*sum(X[i,]*hatPF)-(ki[i]*kbar)/d^2-(stability/2)*((2*kbar*ki[i])/d^2-ki[i]/d-kbar/d+1))
  }
  phi_bar=mean(phi_i)
  var_stab=(4/M^2)*sum((phi_i-phi_bar)^2) ## this is the variance of the stability estimate

  ## then we calculate lower and upper limits of the confidence intervals
  z<-qnorm(1-alpha/2) # this is the standard normal cumulative inverse at a level 1-alpha/2
  upper<-stability+z*sqrt(var_stab) ## the upper bound of the (1-alpha) confidence interval
  lower<-stability-z*sqrt(var_stab) ## the lower bound of the (1-alpha) confidence interval

  return(list("stability"=stability,"variance"=var_stab,"lower"=lower,"upper"=upper))

}


# AirHOLP screening
#
# Performs Air-HOLP variable screening.
#
# AirHOLP (2025) from "https://github.com/Logic314/Air-HOLP"
#
#' @keywords internal
#' @noRd
AirHOLP <- function(X, y, Threshold, r0 = 10, adapt = TRUE,
                    iter = 10, Lambda, U, XU) {
  # Arguments:-
  # X: matrix of features (Matrix)
  # y: response vector (Vector)
  # Threshold: screening threshold (Integer)
  # r0: initial penalties (Vector)
  # adapt: if >= 1 adaptive penalty will be used (Binary)
  # iter: maximum number of iterations for adaptive penalty selection (Integer)
  # Lambda: eigenvalues of XXT, if missing the function will compute it (Vector)
  # U: eigenvectors of XXT, if missing the function will compute it (Matrix)
  # XU: X transpose times U, if missing the function will compute it (Matrix)

  # Output:-
  # index_r: ranking of features by Air-HOLP (Matrix)
  # index_r0: ranking of features by Ridge-HOLP (Matrix)
  # Beta_r: regression coefficients of Air-HOLP (Matrix)
  # Beta_r0: regression coefficients of Ridge-HOLP (Matrix)
  # r: selected penalty parameters by Air-HOLP (Vector)
  # iter_last: number of iterations used in Air-HOLP (Vector)

  n <- dim(X)[1] # sample size
  p <- dim(X)[2] # number of features
  q <- length(r0) # number of penalty parameters
  iter_temp2 <- 0*(1:q) # used for calculating iter_last
  iter_temp1 <- iter_temp2 - 1 # used for calculating iter_last

  # Standardizing X and y:
  X <- X - matrix(rep(colMeans(X),each = n),n,p)
  X <- X/matrix(rep(sqrt(colMeans(X^2)),each = n),n,p)
  y <- (y - mean(y))/sd(y)

  if(adapt){
    # Main computations:
    if(missing(Lambda)|missing(U)){
      XXT <- tcrossprod(X)
      eXXT <- eigen(XXT)
      Lambda <- eXXT$values
      U <- eXXT$vectors
    }
    if(missing(XU)){
      XU <- crossprod(X,U)
    }
    Dn <- diag(Lambda)
    UTy <- crossprod(U,y)
    yUD2UTy <- UTy^2*(Lambda^2)

    # Penalty selection:
    r_max <- 1000*sqrt(n) # maximum penalty
    max.iter <- 30 # maximum number of iterations for Newtons method
    index_r <- matrix(1:(p*q), nrow = p, ncol = q)
    index_r0 <- index_r
    Beta_r <- index_r
    Beta_r0 <- index_r
    r <- r0
    r_temp <- r0
    for (j in 1:iter) {
      for (i in 1:q) {
        # Initial screening:
        Beta_temp <- XU%*%((Lambda+r[i])^(-1)*UTy)
        index_temp <- match(1:p,rank(-abs(Beta_temp), na.last = NA,
                                     ties.method = c("random")))
        Xs <- X[,index_temp[1:Threshold]] # Screened features
        if(j<2) {
          Beta_r0[,i] <- Beta_temp
          index_r0[,i] <- rank(-abs(Beta_temp), na.last = NA,
                               ties.method = c("random"))
        }

        # Estimating the expected response:
        ys <- Xs%*%(solve(crossprod(Xs) +
                            diag(Threshold)*10^-12)%*%crossprod(Xs,y))

        # MSE functions:
        ysUDUTy <- t(crossprod(ys,U)*Lambda)*UTy
        Z <- function(lam) { # The function we minimize
          t((Lambda+c(lam))^-2)%*%yUD2UTy - 2*t((Lambda+c(lam))^-1)%*%ysUDUTy
        }
        Z1 <- function(lam) { # First derivative
          -2*t((Lambda+c(lam))^-3)%*%yUD2UTy + 2*t((Lambda+c(lam))^-2)%*%ysUDUTy
        }
        Z2 <- function(lam) { # Second derivative
          6*t((Lambda+c(lam))^-4)%*%yUD2UTy - 4*t((Lambda+c(lam))^-3)%*%ysUDUTy
        }

        # MSE minimization:
        sol <- newton(Z1, Z2, 0.0001, tol = 0.001, m = max.iter)
        r[i] <- sol
        if(r[i] > r_max) {r[i] <- r_max}
        if(r[i] < 0.0001) {r[i] <- 0.0001}
        if(Z(r_max) < Z(r[i])) {r[i] <- r_max} # Checking boundaries
        if(Z(0.0001) < Z(r[i])) {r[i] <- 0.0001}

        # Feature screening:
        Beta_r[,i] <- XU%*%((Lambda+r[i])^(-1)*UTy)
        index_r[,i] <- rank(-abs(Beta_r[,i]), na.last = NA,
                            ties.method = c("random"))

        # Calculations for the number of iterations:
        if(abs(r[i] - r_temp[i]) < 0.01*r[i]){ # Checking relative error
          iter_temp1[i] <- j
          iter_temp2[i] <- iter_temp2[i] + 1
        }
      }
      if(sum(abs(r - r_temp) < 0.01*r) == q){ # Checking relative error
        break
      }
      r_temp <- r
    }
    iter_last <- iter_temp1 - iter_temp2 + 1 # Number of iterations
    AirHOLP <- list(index_r = index_r, index_r0 = index_r0, Beta_r = Beta_r,
                    Beta_r0 = Beta_r0, r = r, iter_last = iter_last)
  } else{
    if(q < 2) {
      # Feature screening:
      if(missing(Lambda)|missing(U)){
        Beta_r0 <- crossprod(X, solve(tcrossprod(X)+r0*diag(n),y))
      } else{
        UTy <- crossprod(U,y)
        Beta_r0 <- XU%*%((Lambda+r0)^(-1)*UTy)
      }
      index_r0 <- rank(-abs(Beta_r0), na.last = NA, ties.method = c("random"))
      AirHOLP <- list(index_r0 = index_r0, Beta_r0 = Beta_r0)
    } else{
      # Main computations:
      if(missing(Lambda)|missing(U)){
        XXT <- tcrossprod(X)
        eXXT <- eigen(XXT)
        Lambda <- eXXT$values
        U <- eXXT$vectors
      }
      if(missing(XU)){
        XU <- crossprod(X,U)
      }
      Dn <- diag(Lambda)
      UTy <- crossprod(U,y)

      # Feature screening:
      index_r <- matrix(1:(p*q), nrow = p, ncol = q)
      for (i in 1:q) {
        Beta_r0[,i] <- XU%*%((Lambda+r0[i])^(-1)*UTy)
        index_r0[,i] <- rank(-abs(Beta_r0[,i]), na.last = NA,
                             ties.method = c("random"))
      }
      AirHOLP <- list(index_r0 = index_r0, Beta_r0 = Beta_r0)
    }
  }
  return(AirHOLP)
}

# Gram-Schmidt orthogonalization
#
# Performs QR decomposition using Gram-Schmidt.
#
# Gram-Schmidt QR decomposition from "https://stackoverflow.com/questions/15584221/gram-schmidt-with-r"
#
#' @keywords internal
#' @noRd
gram_schmidtR <- function(A) {
  m <- nrow(A)
  n <- ncol(A)
  Q <- matrix(0, nrow = m, ncol = n)
  R <- matrix(0, nrow = n, ncol = n)
  for (j in 1:n) {
    v <- A[ , j, drop = FALSE]
    if (j > 1) {
      for(i in 1:(j-1)) {
        R[i, j] <- t(Q[,i,drop = FALSE]) %*% A[ , j, drop = FALSE]
        v <- v - R[i, j] * Q[ ,i]
      }
    }
    R[j, j] = norm(v, type = "2")
    if (R[j, j] == 0)
      stop("Gram-Schmidt failed because the columns are linearly dependent.")
    Q[ ,j] = v / R[j, j]
  }

  list("Q" = Q, "R" = R)

}


