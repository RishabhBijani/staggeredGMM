# Reduced-space optimal-GMM solve.
#
# The N_2x2 clean comparisons are linear images of only the (G+1)T
# group-by-period means, so Omega = U W U' is rank-deficient with a fixed,
# data-independent column space spanned by the left singular vectors of U.
# Because U, Q_H and Delta do not change across GMM iterations, the SVD and
# the projections
#
#   uQ  = u' Q_H     (m x N_beta)
#   uDe = u' Delta   (m)
#
# are computed ONCE. Each iteration then works entirely in the
# m = (G+1)T dimensional space, and no object of size N_2x2 is ever formed.
#
# With Omega^+ = u M^+ u' and M = diag(d) (v' W v) diag(d), the optimal-GMM
# normal equations are
#
#   Q' Omega^+ Q = uQ' M^+ uQ,     Q' Omega^+ Delta = uQ' M^+ uDe
#
# where M^+ is the truncated inverse retaining eigen-directions above
# tol_rel * lambda_max. This is algebraically identical to optimal GMM on a
# maximal linearly independent subset of the comparisons (Proposition 1).


# One-off factorisation of the fixed incidence factor U.
# @noRd
gmm_reduced_factor <- function(U, Q_H, Delta) {
  sv <- svd(U)
  list(d   = sv$d,
       v   = sv$v,
       uQ  = crossprod(sv$u, Q_H),
       uDe = as.numeric(crossprod(sv$u, Delta)),
       rank = sum(sv$d > max(sv$d) * .Machine$double.eps * max(dim(U))))
}


# Per-iteration solve given the precomputed factor and the current W.
# Returns NULL when no eigen-direction survives truncation.
# @noRd
gmm_reduced_solve <- function(f, W, tol_rel = 1e-8) {
  M <- outer(f$d, f$d) * crossprod(f$v, W %*% f$v)
  M <- (M + t(M)) / 2
  eg <- eigen(M, symmetric = TRUE)
  lam <- eg$values
  lam_max <- max(lam)
  if (!is.finite(lam_max) || lam_max <= 0) return(NULL)
  keep <- lam > tol_rel * lam_max
  if (!any(keep)) return(NULL)

  Vk <- eg$vectors[, keep, drop = FALSE]
  dk <- 1 / lam[keep]
  Mp_uQ  <- Vk %*% (dk * crossprod(Vk, f$uQ))
  Mp_uDe <- Vk %*% (dk * crossprod(Vk, f$uDe))

  list(QtAQ    = crossprod(f$uQ, Mp_uQ),
       QtAD    = as.numeric(crossprod(f$uQ, Mp_uDe)),
       n_kept  = sum(keep),
       n_total = length(lam))
}


# Sandwich variance for the identity-weighted seed, reconstructed in the
# reduced space: (Q'Q)^{-1} Q' Omega Q (Q'Q)^{-1}.
# @noRd
gmm_seed_sandwich <- function(f, W, QtQ_inv) {
  M <- outer(f$d, f$d) * crossprod(f$v, W %*% f$v)
  M <- (M + t(M)) / 2
  QtOmegaQ <- crossprod(f$uQ, M %*% f$uQ)
  QtQ_inv %*% QtOmegaQ %*% QtQ_inv
}
