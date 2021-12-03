#' Performs cross-validation to calculate the average predicted log likelihood for the \code{logreg2ph} method. This function can be used to select the B-spline basis that yields the largest average predicted log likelihood.
#'
#' @param nfolds Specifies the number of cross-validation folds. The default value is \code{5}. Although \code{nfolds} can be as large as the sample size (leave-one-out cross-validation), it is not recommended for large datasets. The smallest value allowable is \code{3}.
#' @param interp Indicator of whether the B-spline coefficients in the testing data should be linearly interpolated from the training data. Defaults to \code{TRUE}.
#' @param seed (For reproducibility in assigning the folds) an integer to specify the random number generator.
#' @param Y_unval Column name with the unvalidated outcome. If \code{Y_unval} is null, the outcome is assumed to be error-free.
#' @param Y_val Column name with the validated outcome.
#' @param X_unval Column name(s) with the unvalidated predictors.  If \code{X_unval} and \code{X_val} are \code{null}, all precictors are assumed to be error-free.
#' @param X_val Column name(s) with the validated predictors. If \code{X_unval} and \code{X_val} are \code{null}, all precictors are assumed to be error-free.
#' @param C (Optional) Column name(s) with additional error-free covariates.
#' @param Validated Column name with the validation indicator. The validation indicator can be defined as \code{Validated = 1} or \code{TRUE} if the subject was validated and \code{Validated = 0} or \code{FALSE} otherwise.
#' @param Bspline Vector of column names containing the B-spline basis functions.
#' @param data A dataframe with one row per subject containing columns: \code{Y_unval}, \code{Y_val}, \code{X_unval}, \code{X_val}, \code{C}, \code{Validated}, and \code{Bspline}.
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param gamma_pred Vector of columns in \code{data} that pertain to the predictors in the outcome error model.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence. Defaults to \code{1E-4}.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm. Defaults to \code{1000}.
#' @return
#' \item{avg_pred_loglike}{Stores the average predicted log likelihood.}
#' \item{pred_loglike}{Stores the predicted log likelihoood in each fold.}
#' \item{converged}{Stores the convergence status of the EM algorithm in each run.}
#' @export
cv_loglik <- function(seed = 1, interp = TRUE, nfolds = 5, Y_unval = NULL, Y_val = NULL, X_unval = NULL, X_val = NULL, C = NULL,
                      Validated = NULL, Bspline = NULL, data, theta_pred = NULL, gamma_pred = NULL,
                      TOL = 1E-4, MAX_ITER = 1000) {
  if (is.null(theta_pred)) { theta_pred <- c(X_val, C) }
  if (is.null(gamma_pred) & !is.null(Y_unval)) { gamma_pred <- c(X_unval, Y_val, X_val, C) }

  set.seed(seed)
  assign_folds <- sample(x = 1:nfolds, size = nrow(data), replace = TRUE)
  status <- rep(TRUE, nfolds)
  msg <- rep("", nfolds)
  ll <- rep(NA, nfolds)
  for (f in 1:nfolds) {
    train <- data[assign_folds != f, ]
    suppressMessages(
      train_fit <- logreg2ph(Y_unval = Y_unval, Y_val = Y_val, X_unval = X_unval, X_val = X_val, C = C,
                             Validated = Validated, Bspline = Bspline, data = train,
                             theta_pred = theta_pred, gamma_pred = gamma_pred,
                             noSE = TRUE, TOL = TOL, MAX_ITER = MAX_ITER)
    )
    status[f] <- train_fit$converged
    msg[f] <- train_fit$converged_msg

    if (train_fit$converged) {
      train_theta <- train_fit$model_coeff$coeff
      train_gamma <- train_fit$outcome_error_coeff$coeff
      train_p <- train_fit$bspline_coeff
      train_x <- data.frame(train[train[, Validated] == 1, X_val])
      train_x <- data.frame(train_x[order(train_x[, 1]), ])
      colnames(train_x) <- X_val
      train_x <- cbind(k = 1:nrow(train_x), train_x)
      train_p <- merge(train_x, train_p)

      test <- data[assign_folds == f, ]

      if (interp) {
        test_x <- data.frame(test[test[, Validated] == 1, X_val])
        test_x <- data.frame(test_x[order(test_x[, 1]), ])
        colnames(test_x) <- X_val
        test_x <- cbind(k_ = 1:nrow(test_x), test_x)
        test_p <- matrix(data = NA, nrow = nrow(test_x), ncol = length(Bspline))

        for (i in 1:nrow(test_x)) {
          x_ <- test_x[i, X_val]
          bf <- suppressWarnings(expr = max(which(train_x[, X_val] <= x_)))
          af <- suppressWarnings(expr = min(which(train_x[, X_val] >= x_)))
          if (bf == -Inf) { bf <- af }
          if (af == Inf) { af <- bf }

          # x values immediately before/after
          x0 <- train_p[bf, X_val]
          x1 <- train_p[af, X_val]

          # B-spline coefficients immediately before/after
          p0 <- train_p[bf, -c(1:(1 + length(X_val)))]
          p1 <- train_p[af, -c(1:(1 + length(X_val)))]

          if (x1 == x0) {
            test_p[i, ] <- unlist(p0)
          } else {
            test_p[i, ] <- unlist((p0 * (x1 - x_) + p1 * (x_ - x0)) / (x1 - x0))
          }
        }

        # Recale columns of test_p to sum to 1
        denom <- colSums(test_p)
        denom[denom == 0] <- 1 # Avoid NaN error due to dividing by 0
        re_test_p <- t(t(test_p) / denom)

        # Construct complete dataset
        cd <- complete_data(Y_unval = Y_unval, Y_val = Y_val, X_unval = X_unval, X_val = X_val, C = C,
                            Validated = Validated, Bspline = Bspline, data = test)
        # Calculate log-likelihood -------------------------------------------
        ll_f <- observed_data_loglik(N = nrow(test), n = sum(test[, Validated]),
                                     Y_unval = Y_unval, Y_val = Y_val, X_unval = X_unval, X_val = X_val, C = C,
                                     Bspline = Bspline, comp_dat_all = cd, theta_pred = theta_pred, gamma_pred = gamma_pred,
                                     theta = train_theta, gamma = train_gamma, p = re_test_p)
      } else {
        # Construct complete dataset
        cd <- complete_data(Y_unval = Y_unval, Y_val = Y_val, X_unval = X_unval, X_val = X_val,
                            try_X = train_x, C = C,
                            Validated = Validated, Bspline = Bspline, data = test)

        # Calculate log-likelihood -------------------------------------------
        ll_f <- observed_data_loglik(N = nrow(test), n = sum(test[, Validated]),
                                     Y_unval = Y_unval, Y_val = Y_val, X_unval = X_unval, X_val = X_val, C = C,
                                     Bspline = Bspline, comp_dat_all = cd, theta_pred = theta_pred, gamma_pred = gamma_pred,
                                     theta = train_theta, gamma = train_gamma, p = train_p[, Bspline])
      }
      ll[f] <- ll_f
    }
  }
  return(list(avg_pred_loglike = mean(ll), pred_loglike = ll, converged = msg))
}

