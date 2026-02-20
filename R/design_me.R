# Design primitives + simulation sampler for the vanishing measurement-error Monte Carlo design.
# The design object is the single source of truth for both simulation and population calculations.
#
# Exports:
#   - make_me_design()
#   - simulate_me_design()
#
# Internal (not exported):
#   - .me_B_matrix()
#   - .me_beta_star_default()
#   - .me_Delta_base_matrix()
#   - .me_proj_range_from_B()
#   - .me_psd_sqrt()

# ---------- internal helpers ----------

# Internal: B matrix for latent factor design
.me_B_matrix <- function(case) {
  if (case == "regular") {
    diag(5)
  } else {
    # k = 3, p = 5
    rbind(
      c(1, 0, 0),
      c(0, 1, 0),
      c(0, 0, 1),
      c(0, 1, 1),
      c(0, 1, -1)
    )
  }
}

# Internal: default beta_star by regime (option C for sparse)
.me_beta_star_default <- function(regime) {
  switch(
    regime,
    dense  = c(1.2, -0.8, 0.6, 1.0, -1.1),
    # your option C: beta* sparse and (in singular/nearly) beta0^rls stays sparse because beta* ∈ Range(Q0)
    sparse = c(1.2, 0.0, 0.9, 0.9, -0.9),
    stop("Unknown regime: ", regime)
  )
}

# Internal: baseline Delta_base = D^{1/2} R(rho) D^{1/2} (PD by construction)
.me_Delta_base_matrix <- function(rho = 1/3) {
  D <- diag(c(1, 1, 1, 3, 4))
  R <- matrix(0, 5, 5)
  diag(R) <- 1

  R[2, 4] <- R[4, 2] <- rho
  R[2, 5] <- R[5, 2] <- rho
  R[3, 4] <- R[4, 3] <- rho
  R[3, 5] <- R[5, 3] <- rho
  R[4, 5] <- R[5, 4] <- rho

  Dsqrt <- diag(sqrt(diag(D)))
  Dsqrt %*% R %*% Dsqrt
}

# Internal: projector onto Range(Q0) when Q0 = B B' and B has full column rank
# P = B (B'B)^{-1} B'
.me_proj_range_from_B <- function(B) {
  BtB_inv <- solve(crossprod(B))  # (k x k)
  B %*% BtB_inv %*% t(B)          # (p x p)
}

# Internal: factor A (possibly rectangular) such that Delta = A A' for PSD Delta
# Uses eigen-decomposition; Delta is symmetric PSD by construction in this design.
.me_psd_sqrt <- function(Delta, tol = 1e-12) {
  e <- eigen(Delta, symmetric = TRUE)
  vals <- e$values
  U <- e$vectors

  scale <- max(abs(vals), 1)
  keep <- vals > tol * scale

  if (!any(keep)) {
    p <- nrow(Delta)
    return(list(A = matrix(0, p, 0), rank = 0L, eigenvalues = numeric(0)))
  }

  U1 <- U[, keep, drop = FALSE]
  v1 <- vals[keep]
  A <- U1 %*% diag(sqrt(v1), nrow = length(v1), ncol = length(v1))

  list(A = A, rank = ncol(A), eigenvalues = v1)
}

# ---------- exported ----------

#' Create vanishing measurement-error design object
#'
#' Constructs deterministic objects for the Monte Carlo design and caches a
#' square-root factor of the measurement-error covariance used in simulation.
#'
#' Design:
#' \deqn{x_i^\star = B f_i,\quad f_i\sim N(0,I_k),\qquad y_i = x_i^{\star\prime}\beta^\star + u_i,\ u_i\sim N(0,\sigma_u^2).}
#' Observed proxy:
#' \deqn{x_i = x_i^\star + \bar\eta_i,\qquad \bar\eta_i \sim N(0,\Delta/m_n).}
#'
#' Regimes (coefficient sparsity):
#' \itemize{
#' \item \code{"dense"}: uses a dense \eqn{\beta^\star}.
#' \item \code{"sparse"}: option C, uses a sparse \eqn{\beta^\star} (and thus sparse \eqn{\beta_0^{rls}}).
#' }
#'
#' Cases (identification geometry):
#' \itemize{
#' \item \code{"regular"}: \eqn{k=5}, \eqn{B=I_5}, \eqn{\Delta=\Delta_{\mathrm{base}}} (PD).
#' \item \code{"singular"}: \eqn{k=3}, \eqn{B=B_s}, \eqn{\Delta=P_{\Range(Q_0)}\Delta_{\mathrm{base}}P_{\Range(Q_0)}} (PSD),
#'       so \eqn{Q_{0n}=Q_0+\Delta/m_n} is singular for all \eqn{n} and \eqn{\Range(Q_{0n})=\Range(Q_0)}.
#' \item \code{"nearly_singular"}: \eqn{k=3}, \eqn{B=B_s}, \eqn{\Delta=\Delta_{\mathrm{base}}} (PD),
#'       so \eqn{Q_{0n}} is nonsingular for each \eqn{n} while \eqn{Q_{0n}\to Q_0}.
#' }
#'
#' @param n Sample size.
#' @param case One of c("regular","singular","nearly_singular").
#' @param regime One of c("dense","sparse").
#' @param m_n Number of replicates. If NULL, defaults to \code{ceiling(sqrt(n))} for all cases.
#' @param beta_star Structural coefficient vector (length 5). If NULL, set by \code{regime}.
#' @param sigma_u Error SD.
#' @param rho Correlation parameter for Delta_base construction.
#' @param tol Numerical tolerance used to compute PSD square roots.
#' @return A list (design object) to be passed to \code{simulate_me_design()} and \code{pop_objects_me()}.
#' @export
make_me_design <- function(n,
                           case = c("regular", "singular", "nearly_singular"),
                           regime = c("dense", "sparse"),
                           m_n = NULL,
                           beta_star = NULL,
                           sigma_u = sqrt(2),
                           rho = 1/3,
                           tol = 1e-12) {
  case <- match.arg(case)
  regime <- match.arg(regime)

  if (is.null(m_n)) m_n <- ceiling(sqrt(n))
  if (is.null(beta_star)) beta_star <- .me_beta_star_default(regime)

  B <- .me_B_matrix(case)
  p <- nrow(B)
  k <- ncol(B)

  Delta_base <- .me_Delta_base_matrix(rho)

  Delta <- if (case == "singular") {
    P_range <- .me_proj_range_from_B(B)
    P_range %*% Delta_base %*% P_range
  } else {
    Delta_base
  }

  Delta_sqrt <- .me_psd_sqrt(Delta, tol = tol)

  list(
    n = as.integer(n),
    case = case,
    regime = regime,
    m_n = m_n,
    beta_star = beta_star,
    sigma_u = sigma_u,
    rho = rho,
    tol = tol,
    B = B,
    p = p,
    k = k,
    Delta_base = Delta_base,
    Delta = Delta,
    Delta_sqrt = Delta_sqrt
  )
}

#' Simulate vanishing measurement-error design
#'
#' Draws one i.i.d. sample \((X,y)\) from the design object created by \code{make_me_design()}.
#'
#' @param design Design object created by \code{make_me_design()}.
#' @return A list with \code{X}, \code{y}, \code{X_star}, \code{Ebar}, \code{F}, \code{u}, and \code{design}.
#' @export
simulate_me_design <- function(design) {
  n <- design$n
  p <- design$p
  k <- design$k
  B <- design$B

  F <- matrix(stats::rnorm(n * k), n, k)
  X_star <- F %*% t(B)

  u <- stats::rnorm(n, sd = design$sigma_u)
  y <- drop(X_star %*% design$beta_star + u)

  # bar{eta}_i ~ N(0, Delta/m_n), generated as Z A'/sqrt(m_n) with Delta = A A'
  if (is.finite(design$m_n) && design$m_n > 0 && design$Delta_sqrt$rank > 0) {
    A <- design$Delta_sqrt$A
    r <- design$Delta_sqrt$rank
    Z <- matrix(stats::rnorm(n * r), n, r)
    Ebar <- Z %*% t(A) / sqrt(design$m_n)
    X <- X_star + Ebar
  } else {
    Ebar <- matrix(0, n, p)
    X <- X_star
  }

  list(
    X = X,
    y = y,
    X_star = X_star,
    Ebar = Ebar,
    F = F,
    u = u,
    design = design
  )
}
