# Helper: MSE + bias^2/var decomposition for vector estimators

#' MSE decomposition for coefficient draws
#'
#' Given draws {beta^(r)}_{r=1}^R (rows of `beta_hat`) and a target beta0,
#' computes
#'   MSE = E||beta_hat - beta0||^2 = ||E beta_hat - beta0||^2 + tr(Var(beta_hat)).
#'
#' @param beta_hat Numeric matrix (R x p), each row one draw of the estimator.
#' @param beta0 Numeric vector (length p), target.
#' @return A list with mse, bias2, var, mean_hat, R_eff.
mse_decomp <- function(beta_hat, beta0) {
  ok <- stats::complete.cases(beta_hat)
  beta_hat <- beta_hat[ok, , drop = FALSE]

  if (nrow(beta_hat) == 0) {
    p <- length(beta0)
    return(list(
      mse = NA_real_, bias2 = NA_real_, var = NA_real_,
      mean_hat = rep(NA_real_, p), R_eff = 0L
    ))
  }

  mean_hat <- colMeans(beta_hat)
  bias2 <- sum((mean_hat - beta0)^2)

  vj <- apply(beta_hat, 2, stats::var)
  var_tot <- sum(vj)

  list(
    mse = bias2 + var_tot,
    bias2 = bias2,
    var = var_tot,
    mean_hat = mean_hat,
    R_eff = nrow(beta_hat)
  )
}

# Internal helper: choose one gamma (or minimize over gamma grid) for an R x p x L array
.mse_decomp_select_gamma <- function(beta_hat_array, beta0, gamma_grid, gamma_select = "min") {
  L <- dim(beta_hat_array)[3]
  mse_vec <- rep(NA_real_, L)

  for (ell in seq_len(L)) {
    mse_vec[ell] <- mse_decomp(beta_hat_array[, , ell, drop = FALSE], beta0)$mse
  }

  if (is.character(gamma_select) && gamma_select == "min") {
    ell_star <- which.min(mse_vec)
  } else {
    ell_star <- which.min(abs(gamma_grid - gamma_select))
  }

  list(
    ell = ell_star,
    gamma = gamma_grid[ell_star],
    mse_vec = mse_vec,
    decomp = mse_decomp(beta_hat_array[, , ell_star, drop = FALSE], beta0)
  )
}
