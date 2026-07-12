
#' @noRd
build_cells <- function(start_yrs, T_max) {
  rows <- vector("list", sum(T_max - start_yrs + 1L))
  id_c <- 0L
  for (g in start_yrs)
    for (tt in g:T_max) {
      id_c <- id_c + 1L
      rows[[id_c]] <- data.frame(g = as.integer(g), t = as.integer(tt),
                                 col_id = id_c)
    }
  do.call(rbind, rows)
}
#' @noRd
enumerate_contrasts_clean <- function(cells, all_cohorts, start_yrs, T_min,
                                      has_nt) {
  buf <- vector("list", 1000000L)
  n   <- 0L

  for (idx in seq_len(nrow(cells))) {
    g  <- cells$g[idx]
    tt <- cells$t[idx]
    fc <- cells$col_id[idx]

    if (g > T_min) {
      if (has_nt) {
        for (s in T_min:(g - 1L)) {
          n <- n + 1L
          buf[[n]] <- list(g = g, t_post = tt, c = 0L, t_pre = s,
                           type = "NT", focal_col = fc,
                           bias_neg_col = NA_integer_,
                           bias_pos_col = NA_integer_)
        }
      }
      for (h in all_cohorts[all_cohorts > tt & all_cohorts != 0L]) {
        for (s in T_min:(g - 1L)) {
          n <- n + 1L
          buf[[n]] <- list(g = g, t_post = tt, c = as.integer(h), t_pre = s,
                           type = "NYT", focal_col = fc,
                           bias_neg_col = NA_integer_,
                           bias_pos_col = NA_integer_)
        }
      }
    }
  }

  buf <- buf[seq_len(n)]
  data.frame(
    g            = vapply(buf, `[[`, integer(1L),   "g"),
    t_post       = vapply(buf, `[[`, integer(1L),   "t_post"),
    c            = vapply(buf, `[[`, integer(1L),   "c"),
    t_pre        = vapply(buf, `[[`, integer(1L),   "t_pre"),
    type         = vapply(buf, `[[`, character(1L), "type"),
    focal_col    = vapply(buf, `[[`, integer(1L),   "focal_col"),
    bias_neg_col = vapply(buf, `[[`, integer(1L),   "bias_neg_col"),
    bias_pos_col = vapply(buf, `[[`, integer(1L),   "bias_pos_col"),
    stringsAsFactors = FALSE
  )
}
#' @noRd
build_Q_H <- function(contrasts_df, N_2x2, N_beta) {
  Q_H <- matrix(0, nrow = N_2x2, ncol = N_beta)
  for (r in seq_len(N_2x2)) {
    Q_H[r, contrasts_df$focal_col[r]] <- 1
    if (contrasts_df$type[r] == "AT") {
      bn <- contrasts_df$bias_neg_col[r]
      bp <- contrasts_df$bias_pos_col[r]
      if (!is.na(bn)) Q_H[r, bn] <- Q_H[r, bn] - 1
      if (!is.na(bp)) Q_H[r, bp] <- Q_H[r, bp] + 1
    }
  }
  Q_H
}