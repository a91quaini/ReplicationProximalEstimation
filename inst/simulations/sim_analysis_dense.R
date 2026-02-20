## Example 1 (DENSE): Vanishing measurement-error Monte Carlo
## - runs population analysis + estimation analysis
## - cases: regular, singular, nearly_singular
## - regime: dense
## - n grid: 10, 50, 100, 1000
## - R = 5000 replications

## -----------------------------
## Global configuration
## -----------------------------
regime <- "dense"

n_grid <- c(10L, 50L, 100L, 1000L)
cases <- c("regular", "singular", "nearly_singular")

R_mc <- 5000L
gamma <- 3/4
rho <- 1/3
sigma_u <- sqrt(2)
nu_fun <- function(n) n^(-3/8)

## output paths (relative to package root; create if missing)
out_dir <- file.path("inst", "simulations", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

rds_out <- file.path(out_dir, paste0("sim_ex_1_", regime, "_sim_out.rds"))
pop_rds <- file.path(out_dir, paste0("sim_ex_1_", regime, "_pop_out.rds"))

## -----------------------------
## 1) Population analysis
## -----------------------------
pop_out <- sim_analysis_pop(
  n_grid = n_grid,
  cases = cases,
  regimes = regime,
  sigma_u = sigma_u,
  rho = rho,
  nu_fun = nu_fun
)

saveRDS(pop_out, pop_rds)

## -----------------------------
## Population plots
## -----------------------------
pop_plot_dir <- file.path(out_dir, paste0("pop_plots_", regime))
if (!dir.exists(pop_plot_dir)) dir.create(pop_plot_dir, recursive = TRUE)

for (cs in cases) {
  plot_sim_analysis_pop(pop_out, case = cs, regime = regime, out_dir = pop_plot_dir)
}

message("Saved population plots to directory: ", pop_plot_dir)
message("Saved population objects to: ", pop_rds)

## -----------------------------
## 2) Estimation analysis
## -----------------------------
sim_out <- sim_analysis_est(
  n_grid = n_grid,
  cases = cases,
  regimes = regime,
  R = R_mc,
  gamma = gamma,
  sigma_u = sigma_u,
  rho = rho,
  nu_fun = nu_fun,
  seed = 123
)

saveRDS(sim_out, rds_out)

est_plot_dir <- file.path(out_dir, paste0("est_plots_", regime))
if (!dir.exists(est_plot_dir)) dir.create(est_plot_dir, recursive = TRUE)

for (cs in cases) {
  plot_sim_analysis_est(sim_out, case = cs, regime = regime, nu_fun = nu_fun, out_dir = est_plot_dir)
}

message("Saved estimation plots to directory: ", est_plot_dir)
message("Saved simulation output to: ", rds_out)

## -----------------------------
## Optional: quick console summary
## -----------------------------
for (cs in cases) {
  summ <- plot_sim_analysis_est(
    sim_out, case = cs, regime = regime, nu_fun = nu_fun, out_dir = est_plot_dir, make_plots = FALSE
  )$summary
  message("\nCase: ", cs, " | Regime: ", regime)
  print(summ)
}

