# ------------------------------------------------------------
# R simulation to verify SURE is (approximately) unbiased
# for spectral estimators where the singular values are
# denoised by 1D fused lasso (total variation) on the spectrum.
# ------------------------------------------------------------

# Packages needed:
# install.packages(c("genlasso", "Matrix"))

library(genlasso)
library(Matrix)

# ---------- Helper: fused-lasso proximal on a vector ----------
# Given a vector y (singular values), compute
#   d* = argmin_d 0.5*||y - d||_2^2 + lambda * ||D d||_1
# where D is the 1D first-difference operator.
# Uses genlasso with X = I.
fused_prox <- function(y, lambda, tol = 1e-8) {
  r <- length(y)
  # 1D first-difference matrix of size (r-1) x r
  # genlasso provides a helper getD1d; we build it explicitly to avoid extra deps
  D <- diff(diag(r), differences = 1)
  fit <- genlasso(y = y, X = diag(r), D = D)
  # Interpolate along path at the requested lambda
  pred <- predict(fit, lambda = lambda)
  as.numeric(pred$fit)
}

# ---------- Helper: count fused plateaus in a vector ----------
# Counts contiguous groups where adjacent entries are (numerically) equal.
count_plateaus <- function(d, tol = 1e-7) {
  r <- length(d)
  if (r == 0) return(0L)
  groups <- 1L
  for (i in 2:r) {
    if (abs(d[i] - d[i - 1]) > tol) groups <- groups + 1L
  }
  groups
}

# ---------- Degrees of freedom for spectral TV shrinker ----------
# Given Y and lambda, compute df via the general spectral formula:
# df = tr(J_d) + (m-n) * sum_i d_i/sigma_i + 2 * sum_{i<j} (sigma_i d_i - sigma_j d_j)/(sigma_i^2 - sigma_j^2)
# with tr(J_d) = number of fused plateaus produced by the TV prox on sigma.

df_spectral_tv <- function(Y, lambda, tol = 1e-10) {
  sv <- svd(Y)
  sigma <- sv$d
  r <- length(sigma)
  # Prox on singular values
  d <- fused_prox(sigma, lambda)
  # (1) Jacobian trace term = number of plateaus
  trJ <- count_plateaus(d)
  # (2) (m-n) term
  mn_term <- (nrow(Y) - ncol(Y)) * sum(d / pmax(sigma, tol))
  # (3) Pairwise spectral-geometry term
  pair_term <- 0
  for (i in 1:(r - 1)) {
    for (j in (i + 1):r) {
      denom <- sigma[i]^2 - sigma[j]^2
      if (abs(denom) < tol) {
        # near-tie: take limiting value using a small jitter
        denom <- sign(denom) * tol
      }
      pair_term <- pair_term + 2 * (sigma[i] * d[i] - sigma[j] * d[j]) / denom
    }
  }
  df <- trJ + mn_term + pair_term
  list(df = df, sigma = sigma, d = d, svd = sv)
}

# ---------- Estimator map (for risk and SURE) ----------
# Build Xhat = U diag(d*) V^T given Y and lambda.
Xhat_spectral_tv <- function(Y, lambda) {
  sv <- svd(Y)
  d <- fused_prox(sv$d, lambda)
  sv$u %*% (d * t(sv$v))  # U diag(d) V^T using recycling: d * t(V)
}

# ---------- Monte Carlo experiment ----------
# Compare SURE with true risk E||Xhat - X0||_F^2 under Gaussian noise.
# SURE = ||Y - Xhat||_F^2 + 2 * tau^2 * df - mn * tau^2
# We average over many trials to verify unbiasedness.

simulate_sure_vs_risk <- function(m = 40, n = 30, rank0 = 10, tau = 1.0,
                                  lambda = 1.0, nrep = 200, seed = 1) {
  set.seed(seed)
  # Build a low-rank ground truth X0 with decaying singular values
  U0 <- qr.Q(qr(matrix(rnorm(m * rank0), m, rank0)))
  V0 <- qr.Q(qr(matrix(rnorm(n * rank0), n, rank0)))
  s0 <- seq(from = 5, to = 0.5, length.out = rank0)
  X0 <- U0 %*% (s0 * t(V0))
  
  mn <- m*n
  risk_vals <- numeric(nrep)
  sure_vals <- numeric(nrep)
  
  for (t in 1:nrep) {
    Z <- matrix(rnorm(mn, sd = tau), m, n)
    Y <- X0 + Z
    
    # Estimator
    Xhat <- Xhat_spectral_tv(Y, lambda)
    
    # True risk (single draw): ||Xhat - X0||_F^2
    risk_vals[t] <- sum((Xhat - X0)^2)
    
    # SURE
    df_out <- df_spectral_tv(Y, lambda)
    df <- df_out$df
    resid <- sum((Y - Xhat)^2)
    sure_vals[t] <- resid + 2*tau^2*df - mn*tau^2
  }
  
  list(
    params = list(m = m, n = n, rank0 = rank0, tau = tau, lambda = lambda, nrep = nrep),
    risk_mean = mean(risk_vals),
    sure_mean = mean(sure_vals),
    risk_sd = sd(risk_vals)/sqrt(nrep),
    sure_sd = sd(sure_vals)/sqrt(nrep),
    mean_diff = mean(sure_vals - risk_vals),
    se_diff = sd(sure_vals - risk_vals)/sqrt(nrep),
    all_risk = risk_vals,
    all_sure = sure_vals
  )
}

# ---------- Optional: run a quick demo over a grid of lambdas ----------
run_grid <- function(lambdas = seq(0, 3, length.out = 10), ...) {
  out <- lapply(lambdas, function(lam) {
    res <- simulate_sure_vs_risk(lambda = lam, ...)
    c(lambda = lam,
      risk_mean = res$risk_mean,
      sure_mean = res$sure_mean,
      diff = res$mean_diff,
      se_diff = res$se_diff)
  })
  do.call(rbind, out)
}

# ---------- Example usage (uncomment to run) ----------
grid <- run_grid(lambdas = seq(0, 3, length.out = 8), m = 40, n = 30, rank0 = 5, tau = 1, nrep = 200)
print(round(grid, 3))
#
single <- simulate_sure_vs_risk(lambda = 1.0, m = 40, n = 30, rank0 = 5, tau = 1, nrep = 500)
print(single[c("risk_mean", "sure_mean", "mean_diff", "se_diff")])
#
# Expectation: mean_diff ~ 0 within a couple of standard errors across lambdas.


# ---------- Grid with standard errors (new) ----------
run_grid_with_se <- function(lambdas = seq(0, 3, length.out = 10), ...) {
  out <- lapply(lambdas, function(lam) {
    res <- simulate_sure_vs_risk(lambda = lam, ...)
    data.frame(
      lambda = lam,
      risk_mean = res$risk_mean,
      sure_mean = res$sure_mean,
      risk_se = res$risk_sd,
      sure_se = res$sure_sd,
      diff = res$mean_diff,
      se_diff = res$se_diff
    )
  })
  do.call(rbind, out)
}

# ---------- Plotting utilities (base R, no ggplot) ----------

# Plot mean true risk vs mean SURE across lambdas with ±2*SE ribbons using base graphics.
plot_grid_results <- function(grid_df) {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  ord <- order(grid_df$lambda)
  x <- grid_df$lambda[ord]
  r_mean <- grid_df$risk_mean[ord]
  s_mean <- grid_df$sure_mean[ord]
  r_se <- grid_df$risk_se[ord]
  s_se <- grid_df$sure_se[ord]
  
  # y-limits covering both series and their ±2SE bands
  ymin <- min(r_mean - 2 * r_se, s_mean - 2 * s_se, na.rm = TRUE)
  ymax <- max(r_mean + 2 * r_se, s_mean + 2 * s_se, na.rm = TRUE)
  
  plot(x, r_mean, type = "n", ylim = c(ymin, ymax),
       xlab = expression(lambda), ylab = "Empirical MSE / SURE (mean across trials)",
       main = "Empirical Risk vs. SURE for spectral TV shrinker")
  
  # Shaded ±2SE bands via polygons
  r_low <- r_mean - 2 * r_se; r_high <- r_mean + 2 * r_se
  s_low <- s_mean - 2 * s_se; s_high <- s_mean + 2 * s_se
  
  polygon(c(x, rev(x)), c(r_low, rev(r_high)), border = NA,
          col = adjustcolor("grey", alpha.f = 0.3))
  polygon(c(x, rev(x)), c(s_low, rev(s_high)), border = NA,
          col = adjustcolor("orchid", alpha.f = 0.3))
  
  # Lines and points
  lines(x, r_mean, lty = 1, lwd = 2)
  points(x, r_mean, pch = 16)
  lines(x, s_mean, lty = 2, lwd = 2)
  points(x, s_mean, pch = 1)
  
  legend("topright",
         legend = c("Empirical Risk (±2 SE)", "SURE (±2 SE)"),
         lty = c(1, 2), pch = c(16, 1), bty = "n")
  mtext("Solid = Monte-Carlo, Dashed = SURE; shaded bands are ±2 SE", line = 0.5, cex = 0.8)
}

# Scatter of per-trial SURE vs true risk for a single lambda.
plot_scatter_sure_vs_risk <- function(sim_res) {
  plot(sim_res$all_risk, sim_res$all_sure,
       xlab = "True risk per trial",
       ylab = "SURE per trial",
       main = "Per-trial SURE vs. true risk",
       pch = 16)
  abline(0, 1, lty = 2)
}

# ---------- Example usage for the plots (uncomment to run) ----------
# set.seed(1)
grid <- run_grid_with_se(lambdas = seq(0, 3, length.out = 12),
                          m = 40, n = 30, rank0 = 5, tau = 1, nrep = 300)
print(round(grid, 3))
plot_grid_results(grid)
#
single <- simulate_sure_vs_risk(lambda = 1.0, m = 40, n = 30, rank0 = 5, tau = 1, nrep = 500)
#plot_scatter_sure_vs_risk(single)

# ---------- One-call wrapper to produce and save the plot (base R only) ----------
# This will:
#  1) run the grid of lambdas,
#  2) draw the risk vs SURE curves with ±2SE bands,
#  3) save a PDF next to your session.
run_and_plot <- function(file = "sure_vs_risk.pdf",
                         lambdas = seq(0, 5, length.out = 20),
                         m = 40, n = 30, rank0 = 5, tau = 1, nrep = 300, seed = 1) {
  set.seed(seed)
  grid <- run_grid_with_se(lambdas = lambdas, m = m, n = n, rank0 = rank0, tau = tau, nrep = nrep)
  # Draw on screen
  plot_grid_results(grid)
  # Save to PDF
  pdf(file, width = 7, height = 5)
  on.exit(dev.off(), add = TRUE)
  plot_grid_results(grid)
  invisible(grid)
}

run_and_plot()  # produces the plot and saves 'sure_vs_risk.pdf'

