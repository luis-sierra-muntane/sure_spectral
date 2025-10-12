sure_svt <- function(lambda, sigma, ...){
# SURE_SVT:   Stein's Unbiased Risk Estimate (SURE) for Singular Value
#             Thresholding (SVT).
#   USAGE:
#       R   =   SURE_SVT(LAMBDA, SIGMA, Y)
#       R   =   SURE_SVT(LAMBDA, SIGMA, S, [M N])    
#       R   =   SURE_SVT(LAMBDA, SIGMA, S, [M N], IS_REAL)
#
#   DESCRIPTION:
#       R   =   SURE_SVT(LAMBDA, SIGMA, Y)
#       returns in R the value of SURE for SVT with threshold LAMBDA >
#       0 evaluated at a matrix Y, which can have real or complex
#       entries. SIGMA > 0 is the standard deviation for the
#       noise.
#       
#       R   =   SURE_SVT(LAMBDA, SIGMA, S, [M N])
#       R   =   SURE_SVT(LAMBDA, SIGMA, S, [M N], IS_REAL)
#       accepts the vector S of observed singular values, and the vector [M N] 
#       containing the number of rows and columns of the observed matrix. The 
#       optional flag IS_REAL indicates whether the observed matrix is real-valued 
#       (if IS_REAL = 1) or not. The default is IS_REAL = 1.
#
#   Adapted from MATLAB implementation in:
#
#       ``Unbiased Risk Estimates for Singular Value Thresholding''
#           E.J.Candes, C.A.Sing-Long, and J.D.Trzasko
  
  args <- list(...)
  if (length(args) == 1) {
    Y <- args[[1]]
    dims <- dim(Y)
    M <- dims[1]; N <- dims[2]
    # Treat as "real" if all imaginary parts are zero
    is_real <- all(Im(Y) == 0)
    s <- svd(Y, nu = 0, nv = 0)$d
  } else if (length(args) == 2) {
    s <- as.numeric(args[[1]])
    M <- as.integer(args[[2]][1])
    N <- as.integer(args[[2]][2])
    is_real <- TRUE
  } else if (length(args) == 3) {
    s <- as.numeric(args[[1]])
    M <- as.integer(args[[2]][1])
    N <- as.integer(args[[2]][2])
    is_real_flag <- args[[3]]
    is_real <- if (is.logical(is_real_flag)) {
      is_real_flag
    } else if (is.numeric(is_real_flag)) {
      is_real_flag != 0
    } else {
      as.logical(is_real_flag)
    }
  } else {
    stop("Usage: sure_svt(lambda, sigma, Y) or sure_svt(lambda, sigma, s, c(M,N), [is_real])")
  }
  
  x <- sure_div_svt(lambda, s, is_real, c(M, N))
  
  risk <- -M * N * sigma^2 + sum(pmin(lambda^2, s^2)) + 2 * sigma^2 * x
  if (!is_real) {
    risk <- risk - M * N * sigma^2
  }
  return(risk)
}

sure_div_svt <- function(lambda, s, is_real, sz) {
  svThreshold <- 1e-8 # threshold to determine whether two singular
                      # values are the same
  M <- as.integer(sz[1]); N <- as.integer(sz[2])
  
  # Robust multiplicity check: build two-column matrix [value, multiplicity]
  if (length(s)==0) stop("Empty singular value vector.")
  z <- if (length(s) >= 2) s[-1] else numeric(0)
  S <- matrix(c(s[1], 1), nrow = 1L, ncol = 2L)
  Is <- 1L
  while (length(z) > 0) {
    idx <- which(abs(z - S[Is, 1]) < svThreshold)
    if (length(idx) == 0) {
      S <- rbind(S, c(z[1], 1))
      z <- z[-1]
      Is <- Is + 1L
    } else {
      z <- z[-idx]
      S[Is, 2] <- S[Is, 2] + length(idx)
    }
  }
  
  # Warnings mirroring the MATLAB prints
  if (any(S[, 1] < svThreshold)) {
    message("   +   [SURE_SVT] Warning: argument might be rank-deficient.")
  }
  if (any(S[, 2] > 1)) {
    message("   +   [SURE_SVT] Warning: argument might have repeated singular values.")
  }
  
  idx_p <- S[, 1] > lambda
  if (isTRUE(is_real)) {
    return(div_svt_real(lambda, S, idx_p, M, N))
  } else {
    return(div_svt_complex(lambda, S, idx_p, M, N))
  }
}

div_svt_real <- function(lambda, S, idx_p, M, N) {
  x <- 0
  if (any(idx_p)) {
    mult <- S[idx_p, 2]
    sval <- S[idx_p, 1]
    x <- x + sum(0.5*mult*(mult + 1))
    x <- x + sum((abs(M - N)*mult + 0.5*mult*(mult - 1))*
                   (pmax(0, sval - lambda)/sval))
  }
  
  K <- nrow(S)
  D <- matrix(0, nrow=K, ncol=K)
  for (Ik in seq_len(K)) {
    denom <- (S[, 1]^2 - S[Ik, 1]^2)
    num <- S[Ik, 2]*S[, 2]*S[, 1]*pmax(0, S[,1] - lambda)
    col <- num/denom
    col[!is.finite(col)|abs(col) > 1e6] <- 0
    D[,Ik] <- col
  }
  
  x <- x + 2*sum(D)
  return(x)
}

div_svt_complex <- function(lambda, S, idx_p, M, N) {
  x <- 0
  if (any(idx_p)) {
    mult <- S[idx_p, 2]
    sval <- S[idx_p, 1]
    x <- x + sum(mult^2)
    x <- x + sum((2 * abs(M - N) + 1 + mult*(mult - 1))*
                   (pmax(0, sval - lambda)/sval))
  }
  
  K <- nrow(S)
  D <- matrix(0, nrow=K, ncol=K)
  for (Ik in seq_len(K)) {
    denom <- (S[, 1]^2 - S[Ik, 1]^2)
    num <- S[Ik, 2]*S[, 2]*S[, 1]*pmax(0, S[, 1] - lambda)
    col <- num/denom
    col[!is.finite(col)|abs(col) > 1e6] <- 0
    D[, Ik] <- col
  }
  
  x <- x + 4 * sum(D)
  return(x)
}