# Population objects for the vanishing measurement-error Monte Carlo design.
#
# Design object (from make_me_design()) must contain:
#   - B, Delta, m_n, beta_star, n, case
#
# Returns population targets:
#   beta0_rls, beta0n_rls, beta0_Delta, beta0n_Delta
# together with Q0, Q0n, delta0, etc.
#
# Conventions:
# - Q0  = B B'
# - Q0n = Q0 + Delta/m_n  (finite m_n), else Q0
# - delta0 = Q0 beta_star
#
# - beta0_rls   := Q0^+ delta0
# - beta0n_rls  := checkQ0n^+ delta0 (eigen hard-threshold at nu)
#
# - beta0n_Delta:
#     * regular, nearly_singular: solve(Q0n, delta0)
#     * singular                : Q0n^+ delta0  (minimum Euclidean norm representative)
#
# - beta0_Delta:
#     * regular        : solve(Q0, delta0)
#     * singular       : Q0^+ delta0  (= beta0_rls)
#     * nearly_singular:
#         - if P0 Delta P0 is nonzero: tie-breaker argmin_{Q0 beta = delta0} beta' Delta beta
#         - else                     : beta0_rls

# -------- helpers (internal) --------

# Moore–Penrose pseudoinverse for symmetric matrices via eigen-decomposition.
.me_pinv_sym <- function(A, tol = 1e-12) {
  e <- eigen(A, symmetric = TRUE)
  vals <- e$values
  U <- e$vectors

  scale <- max(abs(vals), 1)
  keep <- abs(vals) > tol * scale

  if (!any(keep)) {
    return(matrix(0, nrow(A), ncol(A)))
  }

  U1 <- U[, keep, drop = FALSE]
  v1 <- vals[keep]

  U1 %*% (t(U1) / v1)
}

# Orthogonal projector onto Range(Q) for symmetric PSD Q
.me_proj_range <- function(Q, tol = 1e-12) {
  e <- eigen(Q, symmetric = TRUE)
  vals <- e$values
  U <- e$vectors

  scale <- max(abs(vals), 1)
  keep <- vals > tol * scale

  if (!any(keep)) {
    return(matrix(0, nrow(Q), ncol(Q)))
  }

  U1 <- U[, keep, drop = FALSE]
  U1 %*% t(U1)
}

# beta_{0n}^{rls} = checkQ_{0n}^+ delta0 with eigen hard-threshold at nu
.me_beta_rls_finite_n <- function(Q0n, delta0, nu) {
  e <- eigen(Q0n, symmetric = TRUE)
  keep <- e$values >= nu

  if (!any(keep)) {
    return(rep(0, nrow(Q0n)))
  }

  U <- e$vectors[, keep, drop = FALSE]
  vals <- e$values[keep]

  drop(U %*% (drop(crossprod(U, delta0)) / vals))
}

# tie-breaker limit beta_0^Delta = argmin_{Q0 beta = delta0} beta' Delta beta
# Closed form (Delta PD):
#   beta = Delta^{-1} Q0 (Q0 Delta^{-1} Q0)^+ delta0
.me_beta_Delta_limit_tiebreak <- function(Q0, Delta, delta0, tol = 1e-12) {
  invDelta_Q0 <- solve(Delta, Q0)      # Delta^{-1} Q0
  M <- Q0 %*% invDelta_Q0              # Q0 Delta^{-1} Q0
  Mp <- .me_pinv_sym(M, tol = tol)
  drop(invDelta_Q0 %*% (Mp %*% delta0))
}

# -------- exported --------

#' Population objects for the vanishing measurement-error Monte Carlo design
#'
#' Computes:
#' \eqn{\beta_0^{\mathrm{rls}}}, \eqn{\beta_{0n}^{\mathrm{rls}}},
#' \eqn{\beta_0^\Delta}, \eqn{\beta_{0n}^\Delta}.
#'
#' @param design Design object created by \code{make_me_design()}.
#' @param nu Eigenvalue threshold \eqn{\nu_n} for \eqn{\check Q_{0n}}.
#'        Defaults to \code{design$n^(-3/8)}.
#' @param tol Numerical tolerance used in pseudoinverses / rank decisions.
#' @return List with \code{Q0}, \code{Q0n}, \code{Delta}, \code{delta0},
#'         \code{beta0_rls}, \code{beta0n_rls}, \code{beta0_Delta}, \code{beta0n_Delta},
#'         and \code{m_n}, \code{nu}, \code{tol}.
#' @export
pop_objects_me <- function(design,
                           nu = design$n^(-3/8),
                           tol = 1e-12) {

  case <- design$case
  B <- design$B
  Delta <- design$Delta
  m_n <- design$m_n
  beta_star <- design$beta_star

  # symmetric by construction
  Q0  <- tcrossprod(B)  # B B'
  Q0n <- if (is.finite(m_n)) Q0 + Delta / m_n else Q0

  # delta0 = Q0 beta_star (equals delta_{0n} in this design)
  delta0 <- drop(Q0 %*% beta_star)

  # beta_0^{rls} = Q0^+ delta0
  beta0_rls <- drop(.me_pinv_sym(Q0, tol = tol) %*% delta0)

  # beta_{0n}^{rls} = checkQ_{0n}^+ delta0
  beta0n_rls <- .me_beta_rls_finite_n(Q0n, delta0, nu)

  # beta_{0n}^Delta:
  # - regular / nearly_singular: solve
  # - singular: pseudoinverse representative
  beta0n_Delta <- if (case == "singular") {
    drop(.me_pinv_sym(Q0n, tol = tol) %*% delta0)
  } else {
    drop(solve(Q0n, delta0))
  }

  # beta_0^Delta:
  # - regular: solve
  # - singular: pseudoinverse representative (= beta0_rls)
  # - nearly_singular: tie-breaker if noise hits Ker(Q0), else beta0_rls
  beta0_Delta <- if (case == "regular") {
    drop(solve(Q0, delta0))
  } else if (case == "singular") {
    beta0_rls
  } else {
    P_range <- .me_proj_range(Q0, tol = tol)
    P0 <- diag(nrow(Q0)) - P_range
    A <- P0 %*% Delta %*% P0

    if (max(abs(A)) <= tol) {
      beta0_rls
    } else {
      .me_beta_Delta_limit_tiebreak(Q0, Delta, delta0, tol = tol)
    }
  }

  list(
    Q0 = Q0,
    Q0n = Q0n,
    Delta = Delta,
    delta0 = delta0,
    beta0_rls = beta0_rls,
    beta0n_rls = beta0n_rls,
    beta0_Delta = beta0_Delta,
    beta0n_Delta = beta0n_Delta,
    m_n = m_n,
    nu = nu,
    tol = tol
  )
}
