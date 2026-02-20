#' Population analysis under vanishing measurement error (compute-only)
#'
#' For each (case, regime) and each sample size n in `n_grid`, compute the population objects
#' \eqn{\beta_0^{\mathrm{rls}}, \beta_0^{\Delta}, \beta_{0n}^{\mathrm{rls}}, \beta_{0n}^{\Delta}}
#' implied by the measurement-error design. No plots are produced.
#'
#' @param n_grid Integer vector of sample sizes.
#' @param cases Character vector in c("regular","singular","nearly_singular").
#' @param regimes Character vector in c("dense","sparse").
#' @param sigma_u Error SD (passed to `make_me_design()`; irrelevant for population objects).
#' @param rho Correlation parameter for Delta construction.
#' @param nu_fun Function n -> nu_n used for \eqn{\beta_{0n}^{\mathrm{rls}}}. Default n^{-3/8}.
#' @return A list with `by_case` (nested by case then regime) and a long data.frame `df`.
#' @export
sim_analysis_pop <- function(n_grid,
                             cases = c("regular", "singular", "nearly_singular"),
                             regimes = c("dense", "sparse"),
                             sigma_u = sqrt(2),
                             rho = 1/3,
                             nu_fun = function(n) n^(-3/8)) {

  cases <- match.arg(cases, several.ok = TRUE)
  regimes <- match.arg(regimes, several.ok = TRUE)
  n_grid <- as.integer(n_grid)

  out <- vector("list", length(cases)); names(out) <- cases
  df_list <- list(); df_ctr <- 1L

  for (cs in cases) {
    out_cs <- vector("list", length(regimes)); names(out_cs) <- regimes

    for (rg in regimes) {

      # build one design to learn p (and to ensure regime is set by design)
      design0 <- make_me_design(
        n = n_grid[1], case = cs, regime = rg,
        sigma_u = sigma_u, rho = rho
      )
      p <- design0$p

      beta0_rls   <- rep(NA_real_, p)
      beta0_Delta <- rep(NA_real_, p)
      beta0n_rls   <- matrix(NA_real_, length(n_grid), p)
      beta0n_Delta <- matrix(NA_real_, length(n_grid), p)

      for (t in seq_along(n_grid)) {
        n <- n_grid[t]

        design <- make_me_design(
          n = n, case = cs, regime = rg,
          sigma_u = sigma_u, rho = rho
        )

        pop <- pop_objects_me(design, nu = nu_fun(n))

        if (t == 1L) {
          beta0_rls   <- pop$beta0_rls
          beta0_Delta <- pop$beta0_Delta
        }

        beta0n_rls[t, ]   <- pop$beta0n_rls
        beta0n_Delta[t, ] <- pop$beta0n_Delta

        df_list[[df_ctr]] <- data.frame(case = cs, regime = rg, n = n,
                                        object = "beta0n_rls", j = seq_len(p),
                                        value = pop$beta0n_rls)
        df_ctr <- df_ctr + 1L

        df_list[[df_ctr]] <- data.frame(case = cs, regime = rg, n = n,
                                        object = "beta0n_Delta", j = seq_len(p),
                                        value = pop$beta0n_Delta)
        df_ctr <- df_ctr + 1L
      }

      # deviations (L2)
      dev_Delta_n_to_Delta <- rep(NA_real_, length(n_grid))
      dev_rls_n_to_rls <- sqrt(rowSums(
        (beta0n_rls - matrix(beta0_rls, nrow = length(n_grid), ncol = p, byrow = TRUE))^2
      ))
      dev_Delta_n_to_rls <- rep(NA_real_, length(n_grid))

      ok <- is.finite(beta0n_Delta[, 1])
      if (any(ok)) {
        dev_Delta_n_to_Delta[ok] <- sqrt(rowSums(
          (beta0n_Delta[ok, , drop = FALSE] -
             matrix(beta0_Delta, nrow = sum(ok), ncol = p, byrow = TRUE))^2
        ))
        dev_Delta_n_to_rls[ok] <- sqrt(rowSums(
          (beta0n_Delta[ok, , drop = FALSE] -
             matrix(beta0_rls, nrow = sum(ok), ncol = p, byrow = TRUE))^2
        ))
      }

      out_cs[[rg]] <- list(
        n_grid = n_grid,
        beta0_rls = beta0_rls,
        beta0_Delta = beta0_Delta,
        beta0n_rls = beta0n_rls,
        beta0n_Delta = beta0n_Delta,
        dev_beta0n_Delta_to_beta0_Delta = dev_Delta_n_to_Delta,
        dev_beta0n_rls_to_beta0_rls = dev_rls_n_to_rls,
        dev_beta0n_Delta_to_beta0_rls = dev_Delta_n_to_rls
      )

      df_list[[df_ctr]] <- data.frame(case = cs, regime = rg, n = n_grid,
                                      object = "dev_beta0n_Delta_to_beta0_Delta",
                                      j = NA_integer_, value = dev_Delta_n_to_Delta)
      df_ctr <- df_ctr + 1L

      df_list[[df_ctr]] <- data.frame(case = cs, regime = rg, n = n_grid,
                                      object = "dev_beta0n_rls_to_beta0_rls",
                                      j = NA_integer_, value = dev_rls_n_to_rls)
      df_ctr <- df_ctr + 1L

      df_list[[df_ctr]] <- data.frame(case = cs, regime = rg, n = n_grid,
                                      object = "dev_beta0n_Delta_to_beta0_rls",
                                      j = NA_integer_, value = dev_Delta_n_to_rls)
      df_ctr <- df_ctr + 1L
    }

    out[[cs]] <- out_cs
  }

  df <- do.call(rbind, df_list)
  list(by_case = out, df = df)
}

#' Plot population analysis under vanishing measurement error
#'
#' Produces two PNGs for a given (case, regime):
#' (1) grouped barplot comparing beta0^{rls} and beta0^{Delta};
#' (2) L2-deviation curves vs n.
#'
#' Filenames are prefixed by "{case}_{regime}_".
#'
#' @param pop_out Output of \code{sim_analysis_pop()}.
#' @param case One of c("regular","singular","nearly_singular").
#' @param regime One of c("dense","sparse").
#' @param out_dir Directory where PNGs will be written (default ".").
#' @export
plot_sim_analysis_pop <- function(pop_out,
                                  case = c("regular", "singular", "nearly_singular"),
                                  regime = c("dense", "sparse"),
                                  out_dir = ".",
                                  width = 9,
                                  height = 3.2,
                                  dpi = 300) {
  case <- match.arg(case)
  regime <- match.arg(regime)

  res <- pop_out$by_case[[case]][[regime]]

  case_tag <- switch(
    case,
    regular = "regular",
    singular = "singular",
    nearly_singular = "nearlysingular"
  )

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  prefix <- paste0(case_tag, "_", regime)

  ## -----------------------------
  ## 1) Coefficient comparison barplot
  ## -----------------------------
  f1 <- file.path(out_dir, paste0(prefix, "_pop_coeff_comparison_0delta_rls.png"))

  df_bar <- data.frame(
    j = factor(seq_along(res$beta0_rls), levels = seq_along(res$beta0_rls)),
    beta0_rls = as.numeric(res$beta0_rls),
    beta0_Delta = as.numeric(res$beta0_Delta)
  )
  df_bar <- stats::reshape(
    df_bar,
    varying = c("beta0_rls", "beta0_Delta"),
    v.names = "value",
    timevar = "object",
    times = c("beta0_rls", "beta0_Delta"),
    direction = "long"
  )
  df_bar$object <- factor(
    df_bar$object,
    levels = c("beta0_rls", "beta0_Delta"),
    labels = c("rls", "Delta")
  )

  p1 <- ggplot2::ggplot(df_bar, ggplot2::aes(x = j, y = value, fill = object)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.55) +
    ggplot2::scale_fill_manual(
      values = c("blue", "red"),
      name = "Object",
      labels = c(expression(beta[0]^{rls}), expression(beta[0]^Delta))
    ) +
    ggplot2::labs(x = "Coefficient index", y = "Value") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 16),
      legend.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 14),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(filename = f1, plot = p1, width = width, height = height, dpi = dpi)

  ## -----------------------------
  ## 2) L2 deviation curves
  ## -----------------------------
  f2 <- file.path(out_dir, paste0(prefix, "_pop_coeff_l2deviation_0delta_rls.png"))

  df_dev <- data.frame(
    n = as.numeric(res$n_grid),
    dev_Delta_n_to_Delta = as.numeric(res$dev_beta0n_Delta_to_beta0_Delta),
    dev_rls_n_to_rls = as.numeric(res$dev_beta0n_rls_to_beta0_rls),
    dev_Delta_n_to_rls = as.numeric(res$dev_beta0n_Delta_to_beta0_rls)
  )
  df_dev <- stats::reshape(
    df_dev,
    varying = c("dev_Delta_n_to_Delta", "dev_rls_n_to_rls", "dev_Delta_n_to_rls"),
    v.names = "value",
    timevar = "object",
    times = c("Delta_n_to_Delta", "rls_n_to_rls", "Delta_n_to_rls"),
    direction = "long"
  )

  df_dev$color <- factor(
    df_dev$object,
    levels = c("Delta_n_to_Delta", "rls_n_to_rls", "Delta_n_to_rls"),
    labels = c("Delta_n_to_Delta", "rls_n_to_rls", "Delta_n_to_rls")
  )
  df_dev$linetype <- factor(
    df_dev$object,
    levels = c("Delta_n_to_Delta", "rls_n_to_rls", "Delta_n_to_rls"),
    labels = c("dashed", "solid", "dotted")
  )

  p2 <- ggplot2::ggplot(df_dev, ggplot2::aes(x = n, y = value, color = color, linetype = linetype)) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c("red", "blue", "purple"),
      name = "Object",
      labels = c(
        expression(paste("||", beta[0*n]^Delta - beta[0]^Delta, "||")["2"]),
        expression(paste("||", beta[0*n]^{rls} - beta[0]^{rls}, "||")["2"]),
        expression(paste("||", beta[0*n]^Delta - beta[0]^{rls}, "||")["2"])
      )
    ) +
    ggplot2::scale_linetype_identity(guide = "none") +
    ggplot2::labs(x = "n", y = "Value") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 16),
      legend.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 14),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(filename = f2, plot = p2, width = width, height = height, dpi = dpi)

  invisible(res)
}


#' Estimation analysis under vanishing measurement error
#'
#' Simulates data under the measurement-error design and computes, for each replication,
#' the following estimators:
#' (i) ridgeless (RL),
#' (ii) modified ridgeless (MRL),
#' (iii) adaptive ridgeless proximal estimator (RL-AL; ridgeless initial, metric Q_n),
#' (iv) adaptive modified-ridgeless proximal estimator (MRL-AL; modified ridgeless initial,
#'      metric Qbar_n).
#'
#' @param n_grid Integer vector of sample sizes.
#' @param cases Character vector in c("regular","singular","nearly_singular").
#' @param regimes Character vector in c("dense","sparse").
#' @param R Number of Monte Carlo replications per (case, regime, n).
#' @param gamma Exponent for \eqn{\lambda_n = n^{-\gamma}} (default 3/4).
#' @param sigma_u Error SD.
#' @param rho Correlation parameter for Delta construction.
#' @param nu_fun Function n -> nu_n for MRL and for Qbar_n (default n^{-3/8}).
#' @param seed Optional integer seed (set once at function entry).
#' @return Nested list indexed by case, then regime, then n (as character).
#' @export
sim_analysis_est <- function(n_grid,
                             cases = c("regular", "singular", "nearly_singular"),
                             regimes = c("dense", "sparse"),
                             R = 5000L,
                             gamma = 3/4,
                             sigma_u = sqrt(2),
                             rho = 1/3,
                             nu_fun = function(n) n^(-3/8),
                             seed = NULL) {

  cases <- match.arg(cases, several.ok = TRUE)
  regimes <- match.arg(regimes, several.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  n_grid <- as.integer(n_grid)
  R <- as.integer(R)

  out <- vector("list", length(cases)); names(out) <- cases

  for (cs in cases) {
    out_cs <- vector("list", length(regimes)); names(out_cs) <- regimes

    for (rg in regimes) {
      out_rg <- vector("list", length(n_grid))
      names(out_rg) <- as.character(n_grid)

      for (ii in seq_along(n_grid)) {
        n <- n_grid[ii]
        nu <- nu_fun(n)
        lambda <- n^(-gamma)

        design <- make_me_design(
          n = n, case = cs, regime = rg,
          sigma_u = sigma_u, rho = rho
        )
        p <- design$p

        beta_rl    <- matrix(NA_real_, R, p)
        beta_mrl   <- matrix(NA_real_, R, p)
        beta_rlal  <- matrix(NA_real_, R, p)
        beta_mrlal <- matrix(NA_real_, R, p)

        for (r in seq_len(R)) {
          dat <- simulate_me_design(design)
          X <- dat$X
          y <- dat$y

          b_rl  <- est_ridgeless(X, y)
          b_mrl <- est_modified_ridgeless(X, y, nu)

          beta_rl[r, ]  <- b_rl
          beta_mrl[r, ] <- b_mrl

          beta_rlal[r, ] <- tryCatch(
            est_prox_rl_adalasso(X = X, y = y, lambda = lambda, b = b_rl),
            error = function(e) rep(NA_real_, p)
          )

          beta_mrlal[r, ] <- est_prox_mrl_adalasso(
            b = b_mrl, X = X, nu = nu, lambda = lambda
          )
        }

        out_rg[[ii]] <- list(
          design = design,
          nu = nu,
          gamma = gamma,
          lambda = lambda,
          estimates = list(
            rl = beta_rl,
            mrl = beta_mrl,
            rl_adalasso = beta_rlal,
            mrl_adalasso = beta_mrlal
          )
        )
      }

      out_cs[[rg]] <- out_rg
    }

    out[[cs]] <- out_cs
  }

  out
}

#' Plot MSE and bias-variance components for simulated estimators
#'
#' Same as before, but now plots are for a given (case, regime).
#'
#' @param sim_out Output of `sim_analysis_est()`.
#' @param case One of c("regular","singular","nearly_singular").
#' @param regime One of c("dense","sparse").
#' @param nu_fun n -> nu_n used to compute population targets.
#' @param out_dir Output directory for PNG plots.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot resolution.
#' @param make_plots If FALSE, skip plot writing and only return the summary.
#' @return A list with a summary data.frame.
#' @export
plot_sim_analysis_est <- function(sim_out,
                                  case = c("regular", "singular", "nearly_singular"),
                                  regime = c("dense", "sparse"),
                                  nu_fun = function(n) n^(-3/8),
                                  out_dir = ".",
                                  width = 9,
                                  height = 3.2,
                                  dpi = 300,
                                  make_plots = TRUE) {
  case <- match.arg(case)
  regime <- match.arg(regime)

  res_case <- sim_out[[case]][[regime]]
  n_grid <- as.integer(names(res_case))
  o <- order(n_grid)
  n_grid <- n_grid[o]
  res_case <- res_case[o]

  mse_rl <- bias2_rl <- var_rl <- rep(NA_real_, length(n_grid))
  mse_mrl <- bias2_mrl <- var_mrl <- rep(NA_real_, length(n_grid))
  mse_rlal <- bias2_rlal <- var_rlal <- rep(NA_real_, length(n_grid))
  mse_mrlal <- bias2_mrlal <- var_mrlal <- rep(NA_real_, length(n_grid))

  for (i in seq_along(n_grid)) {
    n <- n_grid[i]
    obj <- res_case[[i]]

    design <- obj$design
    pop <- pop_objects_me(design, nu = nu_fun(n))

    # targets: RL -> beta0^Delta ; MRL -> beta0^rls
    target_rl  <- pop$beta0_Delta
    target_mrl <- pop$beta0_rls

    d_rl  <- mse_decomp(obj$estimates$rl,  target_rl)
    d_mrl <- mse_decomp(obj$estimates$mrl, target_mrl)
    mse_rl[i]   <- d_rl$mse
    bias2_rl[i] <- d_rl$bias2
    var_rl[i]   <- d_rl$var
    mse_mrl[i]   <- d_mrl$mse
    bias2_mrl[i] <- d_mrl$bias2
    var_mrl[i]   <- d_mrl$var

    d_rlal  <- mse_decomp(obj$estimates$rl_adalasso,  target_rl)
    d_mrlal <- mse_decomp(obj$estimates$mrl_adalasso, target_mrl)
    mse_rlal[i]   <- d_rlal$mse
    bias2_rlal[i] <- d_rlal$bias2
    var_rlal[i]   <- d_rlal$var
    mse_mrlal[i]   <- d_mrlal$mse
    bias2_mrlal[i] <- d_mrlal$bias2
    var_mrlal[i]   <- d_mrlal$var
  }

  summary_df <- data.frame(
    case = case,
    regime = regime,
    n = n_grid,
    mse_rl = mse_rl, bias2_rl = bias2_rl, var_rl = var_rl,
    mse_mrl = mse_mrl, bias2_mrl = bias2_mrl, var_mrl = var_mrl,
    mse_rlal = mse_rlal, bias2_rlal = bias2_rlal, var_rlal = var_rlal,
    mse_mrlal = mse_mrlal, bias2_mrlal = bias2_mrlal, var_mrlal = var_mrlal
  )

  if (!make_plots) return(list(summary = summary_df))

  case_tag <- switch(
    case,
    regular = "regular",
    singular = "singular",
    nearly_singular = "nearlysingular"
  )

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  prefix <- paste0(case_tag, "_", regime)

  common_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 16),
      legend.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 14),
      panel.grid.minor = ggplot2::element_blank()
    )

  # (1) MSE: RL vs MRL
  df_mse <- data.frame(n = n_grid, RL = mse_rl, MRL = mse_mrl)
  df_mse <- stats::reshape(
    df_mse,
    varying = c("RL", "MRL"),
    v.names = "value",
    timevar = "estimator",
    times = c("RL", "MRL"),
    direction = "long"
  )
  df_mse$estimator <- factor(df_mse$estimator, levels = c("RL", "MRL"))

  p1 <- ggplot2::ggplot(df_mse, ggplot2::aes(x = n, y = value, color = estimator, linetype = estimator)) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(RL = "red", MRL = "blue"), name = "Estimator") +
    ggplot2::scale_linetype_manual(values = c(RL = "dashed", MRL = "solid"), name = "Estimator") +
    ggplot2::labs(x = "n", y = "MSE") +
    common_theme

  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_MSE_comparison_RLMRL.png")),
    plot = p1, width = width, height = height, dpi = dpi
  )

  # (2) Bias^2 and Var: RL vs MRL
  df_comp <- data.frame(
    n = rep(n_grid, times = 4),
    estimator = rep(c("RL", "RL", "MRL", "MRL"), each = length(n_grid)),
    component = rep(c("Bias^2", "Variance", "Bias^2", "Variance"), each = length(n_grid)),
    value = c(bias2_rl, var_rl, bias2_mrl, var_mrl)
  )
  df_comp$estimator <- factor(df_comp$estimator, levels = c("RL", "MRL"))
  df_comp$component <- factor(df_comp$component, levels = c("Bias^2", "Variance"))

  p2 <- ggplot2::ggplot(df_comp, ggplot2::aes(x = n, y = value, color = estimator, linetype = component)) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(RL = "red", MRL = "blue"), name = "Estimator") +
    ggplot2::scale_linetype_manual(values = c("Bias^2" = "dashed", "Variance" = "solid"), name = "Component") +
    ggplot2::labs(x = "n", y = "MSE decomposition") +
    common_theme

  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_BIASVAR_comparison_RLMRL.png")),
    plot = p2, width = width, height = height, dpi = dpi
  )

  # (3) MSE: ARL vs AMRL
  df_mse_adapt <- data.frame(n = n_grid, ARL = mse_rlal, AMRL = mse_mrlal)
  df_mse_adapt <- stats::reshape(
    df_mse_adapt,
    varying = c("ARL", "AMRL"),
    v.names = "value",
    timevar = "estimator",
    times = c("ARL", "AMRL"),
    direction = "long"
  )
  df_mse_adapt$estimator <- factor(df_mse_adapt$estimator, levels = c("ARL", "AMRL"))

  p3 <- ggplot2::ggplot(df_mse_adapt, ggplot2::aes(x = n, y = value, color = estimator, linetype = estimator)) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(ARL = "red", AMRL = "blue"), name = "Estimator") +
    ggplot2::scale_linetype_manual(values = c(ARL = "dashed", AMRL = "solid"), name = "Estimator") +
    ggplot2::labs(x = "n", y = "MSE") +
    common_theme

  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_MSE_comparison_ARLAMRL.png")),
    plot = p3, width = width, height = height, dpi = dpi
  )

  # (4) Bias^2 and Var: ARL vs AMRL
  df_comp_adapt <- data.frame(
    n = rep(n_grid, times = 4),
    estimator = rep(c("ARL", "ARL", "AMRL", "AMRL"), each = length(n_grid)),
    component = rep(c("Bias^2", "Variance", "Bias^2", "Variance"), each = length(n_grid)),
    value = c(bias2_rlal, var_rlal, bias2_mrlal, var_mrlal)
  )
  df_comp_adapt$estimator <- factor(df_comp_adapt$estimator, levels = c("ARL", "AMRL"))
  df_comp_adapt$component <- factor(df_comp_adapt$component, levels = c("Bias^2", "Variance"))

  p4 <- ggplot2::ggplot(df_comp_adapt, ggplot2::aes(x = n, y = value, color = estimator, linetype = component)) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(ARL = "red", AMRL = "blue"), name = "Estimator") +
    ggplot2::scale_linetype_manual(values = c("Bias^2" = "dashed", "Variance" = "solid"), name = "Component") +
    ggplot2::labs(x = "n", y = "MSE decomposition") +
    common_theme

  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_BIASVAR_comparison_ARLAMRL.png")),
    plot = p4, width = width, height = height, dpi = dpi
  )

  list(summary = summary_df)
}
