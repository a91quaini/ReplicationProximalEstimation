#' Ridgeless estimator (sample Moore--Penrose)
#'
#' Implements the ridgeless estimator as in the paper:
#' \deqn{\hat\beta_n^{\mathrm{rl}} = Q_n^+ \delta_n,\quad
#' Q_n = X'X/n,\ \delta_n = X'y/n.}
#' Uses `MASS::ginv()` for the Moore--Penrose pseudoinverse.
#'
#' @param X Numeric matrix (n x p).
#' @param y Numeric vector (length n).
#' @return Numeric vector (length p).
#' @export
est_ridgeless <- function(X, y) {
  n <- nrow(X)
  Qn    <- crossprod(X) / n
  delta <- crossprod(X, y) / n
  drop(MASS::ginv(Qn) %*% delta)
}

#' Modified ridgeless estimator (rank recovery + projection)
#'
#' Implements the modified ridgeless estimator as in the paper:
#' \deqn{\check\beta_n^{ls} = \check Q_n^+ \delta_n,}
#' where \eqn{\check Q_n} is obtained by hard-thresholding the eigenvalues of
#' \eqn{Q_n = X'X/n} at \eqn{\nu}:
#' \deqn{\check Q_n = U \,\mathrm{diag}(\lambda_j 1\{\lambda_j \ge \nu\})\, U'.}
#' Equivalently,
#' \deqn{\check Q_n^+ = U \,\mathrm{diag}(\lambda_j^{-1} 1\{\lambda_j \ge \nu\})\, U'.}
#'
#' @param X Numeric matrix (n x p).
#' @param y Numeric vector (length n).
#' @param nu Eigenvalue threshold \eqn{\nu_n} (nonnegative scalar in the theory).
#' @return Numeric vector (length p).
#' @export
est_modified_ridgeless <- function(X, y, nu) {
  n <- nrow(X)
  Qn    <- crossprod(X) / n
  delta <- crossprod(X, y) / n

  e <- eigen(Qn, symmetric = TRUE)

  keep <- e$values >= nu
  if (!any(keep)) {
    return(rep(0, ncol(X)))
  }

  U    <- e$vectors[, keep, drop = FALSE]
  vals <- e$values[keep]

  # beta = U diag(1/vals) U' delta
  # Compute as: U %*% ((U' delta) / vals)
  drop(U %*% (drop(crossprod(U, delta)) / vals))
}

#' Proximal adaptive lasso based on modified ridgeless, solved via glmnet
#'
#' Implements the adaptive-lasso proximal step with metric \eqn{\overline{\check Q}_n}.
#'
#' @param b Initial estimator (numeric vector of length p). If NULL, it is computed as
#'   \eqn{\check\beta_n^{ls}} from (X,y,nu).
#' @param X Regressor matrix (n x p) used to construct \eqn{Q_n=X'X/n}.
#' @param y Outcome vector (length n). Required only if b is NULL.
#' @param nu Eigenvalue threshold \eqn{\nu_n} used to form \eqn{\check Q_n}.
#' @param lambda Nonnegative tuning parameter \eqn{\lambda_n}.
#' @return Numeric vector of coefficients (length p).
#' @export
est_prox_mrl_adalasso <- function(b = NULL, X, y = NULL, nu, lambda) {
  # If b is not provided, compute the modified ridgeless initial estimator
  if (is.null(b)) {
    if (is.null(y)) stop("est_prox_mrl_adalasso: y must be provided when b is NULL.")
    b <- est_modified_ridgeless(X, y, nu)
  }

  # Prox with lambda=0 returns the input point b
  if (lambda == 0) return(b)

  n <- nrow(X)
  p <- length(b)

  # Build Q_n and its eigendecomposition (Q_n is symmetric by construction)
  Qn <- crossprod(X) / n
  e  <- eigen(Qn, symmetric = TRUE)

  # Retained eigenspace for check Q_n
  # (exclude exact zeros to keep the closed-form Abar strictly PD even if nu=0)
  keep_ev <- (e$values >= nu) & (e$values > 0)

  # Construct Abar = overline{check Q}_n:
  # eigenvalues are {vals_kept} on span(U_kept) and 1 on its orthogonal complement.
  if (!any(keep_ev)) {
    Abar <- diag(p)
  } else {
    U    <- e$vectors[, keep_ev, drop = FALSE]
    vals <- e$values[keep_ev]
    Abar <- diag(p) + U %*% diag(vals - 1, nrow = length(vals)) %*% t(U)
  }

  # Cholesky factor: Abar = R'R (Abar is PD under nu>0; and still PD under keep_ev rule)
  R <- chol(Abar)

  # Adaptive weights: w_j = 1/|b_j| with 1/0 = +Inf -> drop |b_j| = 0
  keep <- abs(b) > 0
  beta_hat <- numeric(p)
  if (!any(keep)) return(beta_hat)

  w_keep <- 1 / abs(b[keep])

  # Weighted lasso reduction: min 0.5||R beta - R b||^2 + lambda * sum w_j |beta_j|
  y_tilde <- drop(R %*% b)
  X_tilde <- R[, keep, drop = FALSE]

  # glmnet objective: (1/(2nobs))||y - Xb||^2 + lambda_glm * sum w|b|
  # Here nobs = nrow(X_tilde) = p, so lambda_glm = lambda / p.
  lambda_glm <- lambda / nrow(X_tilde)

  fit <- glmnet::glmnet(
    x = X_tilde,
    y = y_tilde,
    alpha = 1,
    intercept = FALSE,
    standardize = FALSE,
    penalty.factor = w_keep,
    lambda = lambda_glm
  )

  beta_keep <- drop(glmnet::coef.glmnet(fit, s = lambda_glm))[-1]
  beta_hat[keep] <- beta_keep
  beta_hat
}

#' Proximal adaptive lasso based on ridgeless, solved via glmnet
#'
#' Adaptive-lasso proximal step with ridgeless initial estimator and weight matrix Q_n = X'X/n.
#'
#' Let b denote the initial estimator (typically the ridgeless estimator
#' \eqn{\hat\beta_n^{\mathrm{rl}} = Q_n^+ \delta_n}). Define adaptive weights
#' \eqn{w_j = 1/|b_j|} with the convention \eqn{1/0 = +\infty}. Coordinates with \eqn{|b_j|=0}
#' are removed from the lasso fit (equivalently, assigned infinite weight) and the output
#' sets them to zero.
#'
#' The estimator returned is
#' \deqn{
#' \beta^+
#' \in \arg\min_{\beta\in\mathbb R^p}
#' \left\{
#' \frac12(\beta-b)' Q_n (\beta-b)
#' + \lambda \sum_{j=1}^p \frac{|\beta_j|}{|b_j|}
#' \right\},
#' }
#' where \eqn{Q_n = X'X/n}.
#'
#' Computation uses the reduction: if \eqn{Q_n = R'R} (Cholesky), then the objective equals
#' \eqn{\frac12\|R\beta - Rb\|_2^2 + \lambda \sum_j w_j|\beta_j|}, which is a weighted lasso
#' solved via `glmnet`.
#'
#' @param X Regressor matrix (n x p).
#' @param y Outcome vector (length n).
#' @param lambda Nonnegative tuning parameter.
#' @param b Optional initial estimator. If NULL, it is computed as the ridgeless estimator.
#' @return Numeric vector of coefficients (length p).
#' @export
est_prox_rl_adalasso <- function(X, y, lambda, b = NULL) {
  if (lambda == 0) {
    if (is.null(b)) b <- est_ridgeless(X, y)
    return(b)
  }

  if (is.null(b)) {
    b <- est_ridgeless(X, y)
  }

  n <- nrow(X)
  p <- length(b)

  Qn <- crossprod(X) / n

  # Cholesky factorization (requires Qn SPD; for singular Qn this will fail by construction)
  R <- chol(Qn)

  # Adaptive weights w_j = 1/|b_j| with 1/0 = +Inf -> drop |b_j| = 0
  keep <- abs(b) > 0
  beta_hat <- numeric(p)
  if (!any(keep)) return(beta_hat)

  w_keep <- 1 / abs(b[keep])

  # Weighted lasso reduction: min 0.5||R beta - R b||^2 + lambda * sum w_j |beta_j|
  y_tilde <- drop(R %*% b)
  X_tilde <- R[, keep, drop = FALSE]

  # glmnet objective: (1/(2nobs))||y - Xb||^2 + lambda_glm * sum w|b|
  # Here nobs = nrow(X_tilde) = p, so lambda_glm = lambda / p.
  lambda_glm <- lambda / nrow(X_tilde)

  fit <- glmnet::glmnet(
    x = X_tilde,
    y = y_tilde,
    alpha = 1,
    intercept = FALSE,
    standardize = FALSE,
    penalty.factor = w_keep,
    lambda = lambda_glm
  )

  beta_keep <- drop(glmnet::coef.glmnet(fit, s = lambda_glm))[-1]
  beta_hat[keep] <- beta_keep
  beta_hat
}
