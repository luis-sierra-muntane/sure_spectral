set.seed(1234)

## Experiment parameters
M           <- 200L        # rows
N           <- 500L        # cols
Ns          <- 50L         # Monte Carlo samples
Nl          <- 50L         # number of lambda values
lambda_max  <- 50          # max lambda
SNR         <- c(0.5, 1, 2, 4)

K <- min(M, N)

## Helper to make "jet-like" colors
jet <- grDevices::colorRampPalette(
  c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000")
)

## Frobenius norm
fnorm <- function(A) sqrt(sum(A * A))

## Matrices (means)
## rank(X1) = 200 (full rank)
X1 <- matrix(rnorm(M*N), M, N)
X1 <- X1/fnorm(X1)

## rank(X2) = 0.5*M
X2 <- matrix(rnorm(M*N, 0, 1), nrow = M, ncol = N)
sv  <- svd(X2)
idx <- round(0.5*M)
d2  <- sv$d
if (idx + 1 <= length(d2)) d2[(idx + 1):length(d2)] <- 0
X2  <- sv$u %*% diag(d2, nrow = length(d2), ncol = length(d2)) %*% t(sv$v)
X2  <- X2/fnorm(X2)

## rank(X3) = 0.05*M
X3 <- matrix(rnorm(M*N), M, N)
sv  <- svd(X3, nu = K, nv = K)
idx <- round(0.05 * M)
d3  <- sv$d
if (idx + 1 <= length(d3)) d3[(idx + 1):length(d3)] <- 0
X3  <- sv$u %>% { . %*% diag(d3, length(d3), length(d3)) %*% t(sv$v) }  # no diag(.) piping
X3  <- X3/fnorm(X3)

## X4: sigma_i = sqrt(M)/(1 + exp((i - M/2)/20)), i = 1..M
X4 <- matrix(rnorm(M*N, 0, 1), nrow = M, ncol = N)
sv  <- svd(X4)
d4  <- sqrt(M) / (1 + exp(((1:M) - 0.5 * M) / 20))
X4  <- sv$u %*% diag(d4, nrow = length(d4), ncol = length(d4)) %*% t(sv$v)
X4  <- X4/fnorm(X4)

## Containers
X0      <- list(X1, X2, X3, X4)
lambda  <- vector("list", 16L)
tau_w   <- matrix(0.0, nrow = 4L, ncol = 4L)
MCR     <- array(0.0, dim = c(Nl, 4L, 4L))
MCS     <- array(0.0, dim = c(Nl, 4L, 4L, Ns))
SURE    <- array(0.0, dim = c(Nl, 4L, 4L))
AVG_RANK  <- array(0L, dim = c(Nl, 4L, 4L))
RANKS <- array(0L, dim = c(Nl, 4L, 4L, Ns)) 

## Simulations
for (Ik in 1:4) {
  X <- X0[[Ik]]
  cat("Matrix X^0_", Ik, "\n")
  for (In in 1:4) { # loop over SNR
    tau_w[Ik, In] <- 1/(SNR[In]*sqrt(M*N)) # noise standard-deviation
    lam_vec <- seq(0, lambda_max*tau_w[Ik, In], length.out = Nl) # thresholds
    lambda[[4*(Ik - 1) + In]] <- lam_vec

    for (Il in 1:Nl) {
      lam <- lam_vec[Il]
      
      ## Monte Carlo samples
      cat("Running monte carlo method with lambda =", lam, "\n")
      for (Is in 1:Ns) {
        Y <- X + tau_w[Ik, In]*matrix(rnorm(M*N), M, N)
        
        svy <- svd(Y, nu = K, nv = K)   # skinny SVD
        dth <- pmax(0, svy$d - lam)
        SVT_Y <- svy$u %*% diag(dth, length(dth), length(dth)) %*% t(svy$v)
        
        rcur <- sum(dth > 1e-12)        # store rank
        RANKS[Il, In, Ik, Is] <- rcur
        
        err2 <- sum((SVT_Y - X)^2)
        MCS[Il, In, Ik, Is] <- err2
        MCR[Il, In, Ik] <- MCR[Il, In, Ik] + err2
        AVG_RANK[Il, In, Ik] <- mean(RANKS[Il, In, Ik, 1:Ns])
      }
      
      ## Average Monte Carlo risk
      MCR[Il, In, Ik] <- MCR[Il, In, Ik] / Ns
      
      ## One-shot SURE (single noise draw)
      Y <- X + tau_w[Ik, In]*matrix(rnorm(M*N), M, N)
      SURE[Il, In, Ik] <- sure_svt(lam, tau_w[Ik, In], Y)
    }
  }
}

## Save results
save(SURE, MCR, MCS, SNR, Nl, lambda, tau_w, file = "data_sure_vs_montecarlo.RData")

## Plotting
cmap <- jet(4)
cmap[2] <- rgb(0, 1, 0.5)
cmap[3] <- rgb(1, 0.5, 0)

lfsz <- 20
afsz <- 24
legend_cex <- 0.5

par(mfrow = c(2, 2))
for (In in 1:4) {
  ymax <- 0; xmax <- 0
  for (Ik in 1:4) {
    lam_vec <- lambda[[4 * (Ik - 1) + In]]
    ymax <- max(ymax, MCR[, In, Ik])
    xmax <- max(xmax, lam_vec)
  }
  
  par(pty = "s")
  plot(NA, xlim = c(0, xmax), ylim = c(-0.1 * ymax, 1.05 * ymax),
       xlab = expression(lambda * tau), ylab = "Empirical Risk",
       main = paste0("SNR = ", SNR[In]),
       cex.lab = lfsz/12)
  
  grid(); box()
  
  ## Monte Carlo risk curves
  for (Ik in 1:4) {
    lam_vec <- lambda[[4 * (Ik - 1) + In]]
    lines(lam_vec, MCR[, In, Ik], col = cmap[Ik], lwd = 2)
  }
  
  ## SURE markers (every other lambda)
  pts_idx <- seq(1, Nl, by = 2)
  for (Ik in 1:4) {
    lam_vec <- lambda[[4 * (Ik - 1) + In]]
    points(lam_vec[pts_idx], SURE[pts_idx, In, Ik],
           pch = 3, cex = 1.2, col = cmap[Ik])
  }
  
  legend("bottomright",
         legend = c(expression(X[1]^{0}), expression(X[2]^{0}),
                    expression(X[3]^{0}), expression(X[4]^{0})),
         col = cmap, lwd = 2, pch = 3, pt.cex = 1.1, cex = legend_cex)
}

# Rank–MSE plots (one panel per SNR)
par(mfrow = c(2, 2))
par(mar = c(4,4,2,1))
for (In in 1:4) {
  # axes limits
  xmax <- 0; ymax <- 0
  for (Ik in 1:4) {
    xmax <- max(xmax, AVG_RANK[, In, Ik], na.rm = TRUE)
    ymax <- max(ymax, MCR[,  In, Ik], na.rm = TRUE)
  }
  
  plot(NA, xlim = c(0, xmax * 1.02), ylim = c(0, ymax * 1.05),
       xlab = "rank(SVT(Y))", ylab = "Empirical Risk",
       main = paste0("SNR = ", SNR[In]))
  
  grid(); box()
  
  # draw lines/points for each example, ordered by rank
  for (Ik in 1:4) {
    x <- AVG_RANK[, In, Ik]
    y <- MCR[,  In, Ik]
    ord <- order(x, decreasing = FALSE)  # left->right by rank
    lines(x[ord], y[ord], col = cmap[Ik], lwd = 2)
    points(x[ord], y[ord], col = cmap[Ik], pch = 16, cex = 0.8)
  }
  
  legend("topleft",
         legend = c(expression(X[1]^{0}), expression(X[2]^{0}),
                    expression(X[3]^{0}), expression(X[4]^{0})),
         col = cmap, lwd = 2, pch = 3, pt.cex = 1.1, cex = legend_cex)
}

save(SURE, MCR, MCS, AVG_RANK, RANKS, SNR, Nl, lambda, tau_w,
     file = "data_sure_vs_montecarlo.RData")

## ---------- RANK vs LAMBDA (SVT) ----------
## Settings
avg_over_samples <- FALSE   # TRUE to average ranks over many noise draws
Ns_rank <- 30L              # only used if avg_over_samples = TRUE
set.seed(4321)

K <- min(M, N)
RANK <- array(0.0, dim = c(Nl, 4L, 4L))  # [Il, In, Ik]

for (Ik in 1:4) {
  X <- X0[[Ik]]
  for (In in 1:4) {
    lam_vec <- lambda[[4 * (Ik - 1) + In]]
    tau <- tau_w[Ik, In]
    
    if (!avg_over_samples) {
      ## One fresh noisy draw; SVD once; ranks for all lambda by comparison
      Y <- X + tau * matrix(rnorm(M * N), M, N)
      svy <- svd(Y, nu = 0, nv = 0)  # only singular values needed
      ## rank(SVT(Y; λ)) = # { d_i(Y) > λ }
      RANK[, In, Ik] <- colSums(outer(svy$d, lam_vec, FUN = ">"))
    } else {
      ## Monte-Carlo average rank across Ns_rank draws
      rsum <- rep(0.0, Nl)
      for (Is in 1:Ns_rank) {
        Y <- X + tau * matrix(rnorm(M * N), M, N)
        d <- svd(Y, nu = 0, nv = 0)$d
        rsum <- rsum + colSums(outer(d, lam_vec, FUN = ">"))
      }
      RANK[, In, Ik] <- rsum / Ns_rank
    }
  }
}

## ---------- Plot: rank vs (tau * lambda) ----------
cmap <- grDevices::colorRampPalette(
  c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000")
)(4)
cmap[2] <- rgb(0, 1, 0.5)
cmap[3] <- rgb(1, 0.5, 0)

par(mfrow = c(2, 2))
for (In in 1:4) {
  ## x-limits mirror your error plots convention (tau * lambda on x-axis)
  xmax <- 0
  for (Ik in 1:4) {
    lam_vec <- lambda[[4 * (Ik - 1) + In]]
    xmax <- max(xmax, tau_w[Ik, In] * lam_vec)
  }
  
  par(pty = "s")
  plot(NA,
       xlim = c(0, xmax),
       ylim = c(0, K),
       xlab = expression(lambda * tau),
       ylab = "rank(SVT(Y))",
       main = paste0("SNR = ", SNR[In]),
       cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.0)
  grid(); box()
  
  for (Ik in 1:4) {
    lam_vec <- lambda[[4 * (Ik - 1) + In]]
    lines(tau_w[Ik, In] * lam_vec, RANK[, In, Ik], col = cmap[Ik], lwd = 2)
  }
  
  legend("bottomleft",
         legend = c(expression(X[1]^{0}), expression(X[2]^{0}),
                    expression(X[3]^{0}), expression(X[4]^{0})),
         col = cmap, lwd = 2, pch = 3, pt.cex = 1.1, cex = legend_cex)
}

## -------- Rank vs Degrees of Freedom (SURE sample) --------
## For each example (Ik) and SNR (In), draw one noisy Y,
## then sweep lambda, computing:  rank(SVT(Y; λ))  and  DoF(λ) = div_SVT(Y; λ).

set.seed(2468)

K <- min(M, N)
RANK_vs_DOF <- array(0.0, dim = c(Nl, 4L, 4L))  # y: rank
DOF          <- array(0.0, dim = c(Nl, 4L, 4L))  # x: degrees of freedom

for (Ik in 1:4) {
  X <- X0[[Ik]]
  for (In in 1:4) {
    tau    <- tau_w[Ik, In]
    lamvec <- lambda[[4 * (Ik - 1) + In]]
    
    ## One SURE-style noisy draw (real-valued)
    Y  <- X + tau * matrix(rnorm(M * N), M, N)
    s  <- svd(Y, nu = 0, nv = 0)$d
    
    ## Sweep λ: rank is # { s_i > λ }; DoF is the SVT divergence
    RANK_vs_DOF[, In, Ik] <- vapply(lamvec, function(lam) sum(s > lam), numeric(1))
    DOF[, In, Ik]         <- vapply(lamvec, function(lam) sure_div_svt(lam, s, TRUE, c(M, N)),
                                    numeric(1))
  }
}

## -------- Plot: rank vs DoF (by SNR; colors = examples) --------
op <- options(scipen = 999)   # strongly prefer plain numbers over sci. notation
on.exit(options(op), add = TRUE)

cmap <- grDevices::colorRampPalette(
  c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000")
)(4)
cmap[2] <- rgb(0, 1, 0.5)
cmap[3] <- rgb(1, 0.5, 0)

par(mfrow = c(2, 2), mar = c(4,4,2,1))
for (In in 1:4) {
  xmax <- 0
  for (Ik in 1:4) xmax <- max(xmax, DOF[, In, Ik], na.rm = TRUE)
  
  plot(NA, xlim = c(0, xmax), ylim = c(0, K),
       xlab = "degrees of freedom (divergence)", ylab = "rank(SVT output)",
       main = paste0("SNR = ", SNR[In]),
       cex.axis = 1.0, cex.lab = 1.2, cex.main = 1.2)
  grid(); box()
  
  for (Ik in 1:4) {
    lines(DOF[, In, Ik], RANK_vs_DOF[, In, Ik], col = cmap[Ik], lwd = 2)
    points(DOF[, In, Ik], RANK_vs_DOF[, In, Ik], pch = 16, cex = 0.6, col = cmap[Ik])
  }
  
  legend("bottomright",
         legend = c(expression(X[1]^{0}), expression(X[2]^{0}),
                    expression(X[3]^{0}), expression(X[4]^{0})),
         col = cmap, lwd = 2, pch = 3, pt.cex = 1.1, cex = legend_cex)
}

