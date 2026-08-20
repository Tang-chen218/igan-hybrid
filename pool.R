# 外池 DL 荟萃 + 稿件对标数字。SE：published > CI > constructed（仅 NefIgArd）
EXT_STUDIES <- data.frame(
  author = c(
    "TESTING", "DAPA-CKD IgAN", "NefIgArd", "PROTECT",
    "EMPA-KIDNEY IgAN", "APPLAUSE-IgAN", "ALIGN"
  ),
  n = c(246, 133, 182, 202, 408, 239, 170),
  yi = c(-4.97, -4.70, -5.37, -3.80, -4.12, -6.12, -4.10),
  # Literature SEs (v3.1): 6/7 extractable; NefIgArd still constructed (SD_ext=6/√n)
  sei = c(0.5612, 0.5000, 0.4447, 0.3827, 0.2400, 0.3622, 0.3571),
  se_source = c(
    "ci",          # JAMA Table 2: −4.97 (−6.07, −3.87)
    "published",   # KI: SE = 0.5
    "constructed", # appendix point −5.37; no arm CI/SE
    "ci",          # chronic −3.8 (−4.6, −3.1)
    "published",   # published SE = 0.24
    "ci",          # −6.12 (−6.83, −5.41)
    "ci"           # total −4.1 (−4.8, −3.4)
  ),
  stringsAsFactors = FALSE
)

# Manuscript Table 1 / Supplementary targets (for audit only)
# v3.0: Jiang Elastic g(T) primary; literature SEs; SD_int=6.0; external_data="fixed"
TABLE_CANONICAL <- list(
  seven = list(mu = -4.72, I2 = 80.9, tau2 = 0.607),
  six_excl_applause = list(mu = -4.42, I2 = 52.0, tau2 = 0.160),
  loo = data.frame(
    excluded = c(
      "None (7-study)", "TESTING", "DAPA-CKD IgAN", "NefIgArd", "PROTECT",
      "EMPA-KIDNEY IgAN", "APPLAUSE-IgAN", "ALIGN"
    ),
    mu = c(-4.72, -4.69, -4.72, -4.62, -4.88, -4.84, -4.42, -4.83),
    I2 = c(80.9, 83.8, 84.1, 82.1, 81.2, 80.6, 52.0, 82.9),
    tau2 = c(0.607, 0.677, 0.705, 0.633, 0.623, 0.734, 0.160, 0.720),
    stringsAsFactors = FALSE
  ),
  table2 = list(
    # v3.0: Jiang Elastic + RE-aligned DGP (θ_int~N(μ,τ²)); Methods se_diff 保留
    trad_power = 0.894,
    no_borrow_power = 0.823,
    ep7_power = 0.851,
    ep7_delta = 0.647,
    ep7_doorstop_power = 0.851,
    ep7_doorstop_rate = 0.106,
    ep7_type1_disc0 = 0.028,
    ep7_type1_disc1 = 0.033,
    ep7_2to1_power = 0.830,
    ep7_2to1_noborrow = 0.777,
    ep7_2to1_doorstop = 0.830,
    ep7_3to1_power = 0.787,
    ep7_3to1_noborrow = 0.705,
    ep6_excl_power = 0.896,
    ep6_excl_delta = 0.723,
    ep6_type1_disc0 = 0.028,
    ep6_type1_disc1 = 0.059
  ),
  s4 = data.frame(
    # RE DGP；doorstop 对 Type I 近似冗余（g(T) 已降权）
    discordance = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
    type1_without = c(0.028, 0.032, 0.033, 0.034, 0.031, 0.029, 0.027),
    type1_with = c(0.028, 0.032, 0.033, 0.034, 0.031, 0.029, 0.027),
    doorstop_trigger = c(0.106, 0.165, 0.318, 0.521, 0.718, 0.868, 0.956),
    stringsAsFactors = FALSE
  ),
  sd_sensitivity = data.frame(
    sd = c(5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0),
    trad_power = c(0.970, 0.938, 0.894, 0.841, 0.785, 0.727, 0.672),
    ep_power = c(0.943, 0.901, 0.851, 0.793, 0.735, 0.683, 0.635),
    no_borrow_power = c(0.934, 0.883, 0.823, 0.760, 0.697, 0.637, 0.581),
    stringsAsFactors = FALSE
  )
)

get_script_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }
  ofile <- sys.frames()[[1]]$ofile
  if (!is.null(ofile)) {
    return(dirname(normalizePath(ofile, winslash = "/")))
  }
  getwd()
}

trad_power <- function(n_total, ratio, delta, sd, alpha = 0.05) {
  nc <- floor(n_total / (ratio + 1))
  nt <- n_total - nc
  se <- sd * sqrt(1 / nt + 1 / nc)
  za <- qnorm(1 - alpha / 2)
  1 - pnorm(za - delta / se) + pnorm(-za - delta / se)
}

der_simonian_laird <- function(yi, sei) {
  v <- sei^2
  w <- 1 / v
  mu_fe <- sum(w * yi) / sum(w)
  Q <- sum(w * (yi - mu_fe)^2)
  k <- length(yi)
  df <- k - 1
  c_val <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / c_val)
  wp <- 1 / (v + tau2)
  mu_re <- sum(wp * yi) / sum(wp)
  se_re <- sqrt(1 / sum(wp))
  I2 <- if (Q > 0) max(0, 100 * (Q - df) / Q) else 0

  list(
    mu_fe = mu_fe,
    mu_re = mu_re,
    se_re = se_re,
    tau2 = tau2,
    tau = sqrt(tau2),
    I2 = I2,
    Q = Q,
    df = df,
    ci_lb = mu_re - 1.96 * se_re,
    ci_ub = mu_re + 1.96 * se_re,
    w_re_pct = wp / sum(wp) * 100
  )
}

meta_for_studies <- function(studies) {
  der_simonian_laird(studies$yi, studies$sei)
}

meta_excluding <- function(exclude_authors) {
  idx <- !EXT_STUDIES$author %in% exclude_authors
  list(
    studies = EXT_STUDIES[idx, , drop = FALSE],
    meta = meta_for_studies(EXT_STUDIES[idx, , drop = FALSE])
  )
}

leave_one_out_meta <- function() {
  loo <- lapply(seq_len(nrow(EXT_STUDIES)), function(i) {
    studies <- EXT_STUDIES[-i, , drop = FALSE]
    meta <- meta_for_studies(studies)
    data.frame(
      excluded = EXT_STUDIES$author[i],
      mu_re = meta$mu_re,
      se_re = meta$se_re,
      tau2 = meta$tau2,
      I2 = meta$I2,
      ci_lb = meta$ci_lb,
      ci_ub = meta$ci_ub,
      stringsAsFactors = FALSE
    )
  })
  loo <- do.call(rbind, loo)
  loo$shift <- loo$mu_re - meta_for_studies(EXT_STUDIES)$mu_re
  loo[order(abs(loo$shift), decreasing = TRUE), , drop = FALSE]
}
