
#' @noRd
gmm_reduced_factor <- function(U, Q_H, Delta) {
  sv <- svd(U)
  list(d   = sv$d,
       v   = sv$v,
       uQ  = crossprod(sv$u, Q_H),
       uDe = as.numeric(crossprod(sv$u, Delta)))
}
#' @noRd
gmm_reduced_solve <- function(f, W, tol_rel = 1e-8) {
  M  <- outer(f$d, f$d) * crossprod(f$v, W %*% f$v)
  M  <- (M + t(M)) / 2
  eg <- eigen(M, symmetric = TRUE)
  lam  <- eg$values
  keep <- lam > tol_rel * max(lam)
  if (!any(keep)) return(NULL)
  Vk <- eg$vectors[, keep, drop = FALSE]
  dk <- 1 / lam[keep]
  Mp_uQ  <- Vk %*% (dk * crossprod(Vk, f$uQ))
  Mp_uDe <- Vk %*% (dk * crossprod(Vk, f$uDe))
  list(QtAQ = crossprod(f$uQ, Mp_uQ),
       QtAD = as.numeric(crossprod(f$uQ, Mp_uDe)))
}