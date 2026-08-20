# Monte Carlo：Jiang elastic prior 混合设计
if (!exists("EXT_STUDIES", inherits = TRUE)) {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  sd <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
  source(file.path(sd, "pool.R"), local = FALSE)
}

trial_group_sizes <- function(n_total, ratio) {
  nc <- floor(n_total / (ratio + 1))
  list(n_control = nc, n_treatment = n_total - nc)
}

# T = z^2；g(T)=1/(1+exp(a+b*log T))；在 Tq0 / Tq1 上解 a,b 使 g=C1 / C2
calibrate_elastic_gT <- function(se_int_ctrl, se_ext_pooled2, tau2,
                                 delta_c = 1.5, C1 = 0.99, C2 = 0.01,
                                 q0 = 0.5, q1 = 0.5, R = 5000L, seed = 1L) {
  set.seed(seed)
  se_diff <- sqrt(se_int_ctrl^2 + se_ext_pooled2 + tau2)
  T0 <- rnorm(R, 0, 1)^2
  T_plus <- rnorm(R, delta_c / se_diff, 1)^2
  T_minus <- rnorm(R, -delta_c / se_diff, 1)^2
  T_q0 <- as.numeric(quantile(T0, q0))
  T_q1 <- min(as.numeric(quantile(T_plus, q1)), as.numeric(quantile(T_minus, q1)))
  if (!(T_q1 > T_q0 && T_q0 > 0)) T_q1 <- max(T_q0 * 4, T_q0 + 1e-6)
  logit <- function(p) log(p / (1 - p))
  b <- (logit(C2) - logit(C1)) / (log(T_q0) - log(T_q1))
  a <- -logit(C1) - b * log(T_q0)
  list(a = a, b = b, T_q0 = T_q0, T_q1 = T_q1, delta_c = delta_c, se_diff = se_diff)
}

g_elastic_T <- function(T, a, b, eps = 1e-12) {
  1 / (1 + exp(a + b * log(pmax(T, eps))))
}

.hybrid_iterate <- function(y_ext, se_ext, tau2, y_int_ctrl, y_int_tx,
                            se_int_ctrl, se_int_tx, robust_mix, doorstop,
                            doorstop_threshold, elastic_cal) {
  w_ext <- 1 / (se_ext^2 + tau2)
  y_ext_pooled <- sum(w_ext * y_ext) / sum(w_ext)
  se_ext_pooled2 <- 1 / sum(w_ext)
  z_signed <- (y_int_ctrl - y_ext_pooled) / sqrt(se_int_ctrl^2 + se_ext_pooled2 + tau2)
  delta_ep <- g_elastic_T(z_signed^2, elastic_cal$a, elastic_cal$b) * (1 - robust_mix)
  doorstop_hit <- doorstop && abs(z_signed) > doorstop_threshold
  d_use <- if (doorstop_hit) 0 else delta_ep
  if (d_use < 1e-8) {
    th <- y_int_tx - y_int_ctrl
    se <- sqrt(se_int_tx^2 + se_int_ctrl^2)
  } else {
    prec_int <- 1 / se_int_ctrl^2
    prec_ext <- d_use / (se_ext_pooled2 + tau2)
    mu_ctrl <- (prec_int * y_int_ctrl + prec_ext * y_ext_pooled) / (prec_int + prec_ext)
    th <- y_int_tx - mu_ctrl
    se <- sqrt(se_int_tx^2 + 1 / (prec_int + prec_ext))
  }
  list(
    reject = pnorm(th / max(se, 1e-10)) > 0.975,
    delta_post = d_use,
    doorstop_hit = doorstop_hit
  )
}

# 主分析：fixed 外池 + θ_int ~ N(μ_RE + discordance, τ²)
simulate_hybrid <- function(
    delta_true = 2.0,
    discordance = 0.0,
    studies = EXT_STUDIES,
    exclude_authors = character(),
    tau2 = NULL,
    mu_re = NULL,
    n_total = 300,
    ratio = 1.0,
    sd_int = 6.0,
    n_sim = 10000,
    seed = 42,
    robust_mix = 0.10,
    doorstop = FALSE,
    doorstop_threshold = 1.5,
    delta_c = 1.5,
    external_data = c("fixed", "resimulate")
) {
  external_data <- match.arg(external_data)
  if (length(exclude_authors)) {
    studies <- studies[!studies$author %in% exclude_authors, , drop = FALSE]
  }
  if (nrow(studies) < 2) stop("External pool must contain at least 2 studies.")

  meta <- meta_for_studies(studies)
  tau2_use <- max(if (is.null(tau2)) meta$tau2 else tau2, 1e-10)
  mu_re_use <- if (is.null(mu_re)) meta$mu_re else mu_re

  grp <- trial_group_sizes(n_total, ratio)
  nc <- grp$n_control
  nt <- grp$n_treatment
  se_int_ctrl <- sd_int / sqrt(nc)
  se_int_tx <- sd_int / sqrt(nt)
  se_ext <- studies$sei
  k <- nrow(studies)

  w0 <- 1 / (se_ext^2 + tau2_use)
  elastic_cal <- calibrate_elastic_gT(se_int_ctrl, 1 / sum(w0), tau2_use, delta_c = delta_c)

  set.seed(seed)
  reject <- logical(n_sim)
  delta_post <- numeric(n_sim)
  doorstop_hits <- logical(n_sim)

  for (i in seq_len(n_sim)) {
    y_ext <- if (external_data == "fixed") {
      studies$yi
    } else {
      rnorm(k, rnorm(k, mu_re_use, sqrt(tau2_use)), se_ext)
    }
    mu_int_true <- rnorm(1L, mu_re_use + discordance, sqrt(tau2_use))
    out <- .hybrid_iterate(
      y_ext, se_ext, tau2_use,
      mean(rnorm(nc, mu_int_true, sd_int)),
      mean(rnorm(nt, mu_int_true + delta_true, sd_int)),
      se_int_ctrl, se_int_tx, robust_mix, doorstop, doorstop_threshold, elastic_cal
    )
    reject[i] <- out$reject
    delta_post[i] <- out$delta_post
    doorstop_hits[i] <- out$doorstop_hit
  }

  list(
    rate = mean(reject),
    delta_mean = mean(delta_post),
    doorstop_rate = mean(doorstop_hits),
    n_sim = n_sim,
    seed = seed
  )
}
