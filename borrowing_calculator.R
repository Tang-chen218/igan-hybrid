#!/usr/bin/env Rscript
# ============================================================================
# IgAN-Hybrid — interactive calculator for hybrid IgAN trial designs
# Elastic-prior borrowing from an external control pool (v3.0 engine).
#
# Usage:
#   Rscript borrowing_calculator.R                          # default: 7-study, 1:1, n=300
#   Rscript borrowing_calculator.R --ratio 2to1 --n-sim 2000
#   Rscript borrowing_calculator.R --exclude APPLAUSE-IgAN --ratio 1to1
#   Rscript borrowing_calculator.R --n-total 320 --sd-int 6.5 --delta-tx 2.0
#   Rscript borrowing_calculator.R --method none            # analytic, no MC
#   Rscript borrowing_calculator.R --interactive            # menu-driven
#   Rscript borrowing_calculator.R --export results.csv
#
# Options:
#   --ratio 1to1|2to1|3to1     allocation ratio (default 1to1)
#   --n-total <int>            hybrid total sample size (default 300)
#   --n-conventional <int>     conventional 1:1 total n for comparison (default 370)
#   --sd-int <num>             planning SD, mL/min/1.73m2/yr (default 6.0)
#   --delta-tx <num>           treatment effect (default 2.0)
#   --discordance <num>        control-mean shift for type I (default 1.0)
#   --exclude <author>[,<author>...]  studies to drop from the pool
#   --method elastic|none      borrowing method (default elastic)
#   --n-sim <int>              MC iterations (default 2000)
#   --seed <int>               RNG seed (default 42)
#   --export <file>            write main result row to CSV
#   --interactive              menu-driven input
#   --no-compare               skip the 3-ratio quick comparison
#   --help
#
# In RStudio: source("borrowing_calculator.R") then call calculate_borrowing().
# Requires only base R plus the pool/simulate engines in this directory.
# ============================================================================

# ---- 0. locate and source engines ----
.sd <- {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) {
    dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
  } else {
    of <- sys.frames()[[1]]$ofile
    if (!is.null(of)) dirname(normalizePath(of, winslash = "/")) else getwd()
  }
}
source(file.path(.sd, "pool.R"), local = FALSE)
source(file.path(.sd, "simulate.R"), local = FALSE)

# ---- 1. helpers ----
parse_ratio <- function(ratio) {
  r <- switch(tolower(ratio),
    "1to1" = 1, "1:1" = 1, "2to1" = 2, "2:1" = 2,
    "3to1" = 3, "3:1" = 3,
    stop("Unknown ratio: ", ratio, " (use 1to1, 2to1, or 3to1)"))
  r
}

ratio_label <- function(r) paste0(r, ":1")

# 当前池内 |yi − μ_RE| 最大值；换池后 discordance 应对齐此量级，而非写死 1.0
pool_max_abs_departure <- function(st, mu_re) {
  if (nrow(st) < 1L) return(NA_real_)
  max(abs(st$yi - mu_re))
}

# Power vs SD_internal（对齐 Figure S2 三线）
power_vs_sd <- function(ratio = "1to1", n_total = 300, n_conventional = 370,
                        delta_tx = 2.0, exclude = character(), method = "elastic",
                        n_sim = 2000, seed = 42, sd_grid = seq(5, 8, by = 0.5)) {
  r <- parse_ratio(ratio)
  st <- EXT_STUDIES[!EXT_STUDIES$author %in% exclude, , drop = FALSE]
  if (nrow(st) < 2L) stop("External pool must retain at least 2 studies.")
  meta <- meta_for_studies(st)
  n_sim <- as.integer(n_sim)
  hybrid <- vapply(sd_grid, function(sd) {
    if (method == "elastic") {
      simulate_hybrid(delta_true = delta_tx, discordance = 0,
                      exclude_authors = exclude, tau2 = meta$tau2, mu_re = meta$mu_re,
                      n_total = n_total, ratio = r, sd_int = sd,
                      n_sim = n_sim, seed = seed, external_data = "fixed",
                      delta_c = 1.5)$rate
    } else {
      trad_power(n_total, r, delta_tx, sd)
    }
  }, numeric(1))
  data.frame(
    sd_int = sd_grid,
    hybrid = hybrid,
    no_borrow = vapply(sd_grid, function(sd) trad_power(n_total, r, delta_tx, sd), numeric(1)),
    conventional = vapply(sd_grid, function(sd) trad_power(n_conventional, 1, delta_tx, sd), numeric(1)),
    stringsAsFactors = FALSE
  )
}

# ---- 2. core function ----
calculate_borrowing <- function(ratio = "1to1", n_total = 300, n_conventional = 370,
                                sd_int = 6.0, delta_tx = 2.0, discordance = 1.0,
                                exclude = character(),
                                method = "elastic", n_sim = 2000, seed = 42,
                                compare = TRUE) {
  r <- parse_ratio(ratio)
  n_conventional <- as.integer(n_conventional)
  if (is.na(n_conventional) || n_conventional < 20) {
    stop("n_conventional must be an integer >= 20.")
  }
  discordance <- as.numeric(discordance)
  if (is.na(discordance) || discordance < 0) {
    stop("discordance must be a number >= 0.")
  }
  st <- EXT_STUDIES[!EXT_STUDIES$author %in% exclude, , drop = FALSE]
  if (nrow(st) < 2) stop("External pool must retain at least 2 studies.")
  meta <- meta_for_studies(st)
  max_dep <- pool_max_abs_departure(st, meta$mu_re)

  run <- function(delta, disc) {
    simulate_hybrid(delta_true = delta, discordance = disc,
                    exclude_authors = exclude, tau2 = meta$tau2, mu_re = meta$mu_re,
                    n_total = n_total, ratio = r, sd_int = sd_int,
                    n_sim = n_sim, seed = seed, external_data = "fixed",
                    delta_c = 1.5)
  }

  if (method == "elastic") {
    pw  <- run(delta_tx, 0)
    t0  <- run(0, 0)
    t1  <- run(0, discordance)
    power   <- pw$rate
    type1_0 <- t0$rate
    type1_1 <- t1$rate
    delta_mean <- pw$delta_mean
  } else if (method == "none") {
    power   <- trad_power(n_total, r, delta_tx, sd_int)
    type1_0 <- 0.05
    type1_1 <- NA
    delta_mean <- NA
  } else {
    stop("Unknown method: ", method, " (use elastic or none)")
  }

  grp <- trial_group_sizes(n_total, r)
  placebo_n <- grp$n_control
  # 传统对照：固定 1:1，安慰剂 n = floor(n_conventional/2)
  trad_placebo <- floor(n_conventional / 2)
  no_borrow_power <- trad_power(n_total, r, delta_tx, sd_int)
  trad_pw <- trad_power(n_conventional, 1, delta_tx, sd_int)

  out <- list(
    pool_n_studies = nrow(st), pool_n = sum(st$n),
    meta_mu = meta$mu_re, meta_tau2 = meta$tau2, meta_I2 = meta$I2,
    ratio_label = ratio_label(r), n_total = n_total, n_treatment = grp$n_treatment,
    n_placebo = placebo_n, n_conventional = n_conventional,
    n_conventional_placebo = trad_placebo,
    sd_int = sd_int, delta_tx = delta_tx, discordance = discordance,
    pool_max_departure = max_dep, method = method,
    power = power, delta_mean = delta_mean,
    type1_disc0 = type1_0, type1_disc1 = type1_1,
    no_borrow_power = no_borrow_power, trad_power = trad_pw,
    borrowing_gain = power - no_borrow_power,
    placebo_reduction = 1 - placebo_n / trad_placebo,
    power_vs_trad = power - trad_pw
  )

  if (compare && method != "none") {
    cmp <- lapply(c(1, 2, 3), function(rr) {
      p <- simulate_hybrid(delta_true = delta_tx, discordance = 0,
                           exclude_authors = exclude, tau2 = meta$tau2,
                           mu_re = meta$mu_re, n_total = n_total, ratio = rr,
                           sd_int = sd_int, n_sim = n_sim, seed = seed,
                           external_data = "fixed", delta_c = 1.5)$rate
      g <- trial_group_sizes(n_total, rr)
      data.frame(ratio = paste0(rr, ":1"), power = p, placebo_n = g$n_control,
                 placebo_pct = 1 - g$n_control / trad_placebo,
                 power_vs_trad = p - trad_pw)
    })
    out$comparison <- do.call(rbind, cmp)
  } else {
    out$comparison <- NULL
  }
  out$sd_curve <- power_vs_sd(
    ratio = ratio, n_total = n_total, n_conventional = n_conventional,
    delta_tx = delta_tx, exclude = exclude, method = method,
    n_sim = n_sim, seed = seed
  )
  out
}

# ---- 3. printing ----
pct1 <- function(x) sprintf("%.1f%%", 100 * x)
fmt <- function(x) ifelse(is.na(x), "—", pct1(x))

print_results <- function(res) {
  cat("\n===== RESULTS =====\n")
  cat(sprintf("External Evidence Pool: %d studies (N=%d)\n", res$pool_n_studies, res$pool_n))
  cat(sprintf("RE Meta-Analysis: μ = %.2f, τ² = %.3f, I² = %.1f%%\n",
              res$meta_mu, res$meta_tau2, res$meta_I2))
  cat(sprintf("\nDesign: %s, n=%d (%d Tx : %d Pl)\n", res$ratio_label, res$n_total,
              res$n_treatment, res$n_placebo))
  cat(sprintf("Method: %s | SD_int = %.1f, Δ = %.1f\n", res$method, res$sd_int, res$delta_tx))
  cat("\n  Hybrid Power:          ", fmt(res$power),
      if (res$method == "none") " (analytic)" else " (MC)", sep = "")
  cat("\n  No-Borrowing Power:    ", pct1(res$no_borrow_power), " (analytic)", sep = "")
  cat(sprintf("\n  Traditional 1:1 (n=%d):", res$n_conventional), pct1(res$trad_power),
      " (analytic)", sep = "")
  cat("\n  Borrowing Gain:        +", sprintf("%.1f pp", 100 * res$borrowing_gain), sep = "")
  cat("\n  Placebo Reduction:     ",
      sprintf("%.0f%% (%d vs %d)", 100 * res$placebo_reduction, res$n_placebo,
              res$n_conventional_placebo), sep = "")
  cat("\n  Type I (means equal):  ", fmt(res$type1_disc0), sep = "")
  cat(sprintf("\n  Type I (diff %.1f):     ", res$discordance), fmt(res$type1_disc1), sep = "")
  if (!is.na(res$pool_max_departure)) {
    cat(sprintf("\n  (Pool max |yi−μ|: %.2f; manuscript used 1.0 for the full completed-trials pool)",
                res$pool_max_departure))
  }
  cat("\n  Power vs Traditional:  ", sprintf("%+.1f pp", 100 * res$power_vs_trad), sep = "")
  cat("\n")

  if (!is.null(res$comparison)) {
    cat("\n===== QUICK COMPARISON (", res$pool_n_studies, "-study pool, SD=",
        res$sd_int, ", Δ=", res$delta_tx, ") =====", sep = "")
    cat("\n  Design  Power    Placebo    vs Trad")
    cat("\n  "); cat(paste(rep("-", 44), collapse = "")); cat("\n")
    for (i in seq_len(nrow(res$comparison))) {
      row <- res$comparison[i, ]
      cat(sprintf("  %-6s  %-7s  %-9s  %+.1f pp\n",
                  row$ratio, pct1(row$power),
                  sprintf("%d (%s)", row$placebo_n,
                          sprintf("%.0f%%", 100 * row$placebo_pct)),
                  100 * row$power_vs_trad))
    }
  }
  cat("\n")
  invisible(res)
}

# ---- 4. CLI parsing (base R, no optparse) ----
parse_cli <- function(args) {
  opts <- list(ratio = "1to1", n_total = 300, n_conventional = 370, sd_int = 6.0,
               delta_tx = 2.0, discordance = 1.0, exclude = character(),
               method = "elastic", n_sim = 2000,
               seed = 42, export = NULL, interactive = FALSE, compare = TRUE,
               help = FALSE)
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    val <- function() { i <<- i + 1; if (i > length(args)) stop("Missing value for ", a); args[i] }
    switch(a,
      "--ratio"    = { opts$ratio <- val() },
      "--n-total"  = { opts$n_total <- as.integer(val()) },
      "--n-conventional" = { opts$n_conventional <- as.integer(val()) },
      "--sd-int"   = { opts$sd_int <- as.numeric(val()) },
      "--delta-tx" = { opts$delta_tx <- as.numeric(val()) },
      "--discordance" = { opts$discordance <- as.numeric(val()) },
      "--exclude"  = { opts$exclude <- c(opts$exclude, strsplit(val(), ",")[[1]]) },
      "--method"   = { opts$method <- val() },
      "--n-sim"    = { opts$n_sim <- as.integer(val()) },
      "--seed"     = { opts$seed <- as.integer(val()) },
      "--export"   = { opts$export <- val() },
      "--interactive" = { opts$interactive <- TRUE },
      "--no-compare"  = { opts$compare <- FALSE },
      "--help"     = { opts$help <- TRUE },
      stop("Unknown option: ", a)
    )
    i <- i + 1
  }
  opts
}

# ---- 5. interactive menu ----
interactive_menu <- function() {
  cat("\n========================================\n")
  cat("  IgAN-Hybrid Trial Design Calculator\n")
  cat("========================================\n\n")
  cat("External Evidence Pool (", nrow(EXT_STUDIES), " studies):\n", sep = "")
  for (j in seq_len(nrow(EXT_STUDIES))) {
    s <- EXT_STUDIES[j, ]
    cat(sprintf("  [%d] %-22s yi=%.2f, SE=%.2f\n", j, s$author, s$yi, s$sei))
  }
  cat("\nSelect studies to include (comma-separated; blank = all):\n> ")
  sel <- readline()
  exclude <- character()
  if (nzchar(trimws(sel))) {
    keep <- as.integer(strsplit(sel, "[, ]+")[[1]])
    exclude <- EXT_STUDIES$author[!seq_len(nrow(EXT_STUDIES)) %in% keep]
  }
  pick1 <- function(prompt, n) {
    v <- as.integer(readline(prompt = prompt))
    if (is.na(v) || v < 1 || v > n) 1 else v
  }
  ratio <- c("1to1", "2to1", "3to1")[pick1("Allocation ratio [1] 1:1  [2] 2:1  [3] 3:1  > ", 3)]
  method <- c("elastic", "none")[pick1("Borrowing method [1] Elastic  [2] No Borrowing  > ", 2)]
  n_total <- as.integer(readline(prompt = "Hybrid total sample size [300]: "))
  if (is.na(n_total)) n_total <- 300
  n_conventional <- as.integer(readline(prompt = "Conventional 1:1 total n [370]: "))
  if (is.na(n_conventional)) n_conventional <- 370
  sd_int <- as.numeric(readline(prompt = "Planning SD [6.0]: "))
  if (is.na(sd_int)) sd_int <- 6.0
  delta_tx <- as.numeric(readline(prompt = "Treatment effect Δ [2.0]: "))
  if (is.na(delta_tx)) delta_tx <- 2.0
  discordance <- as.numeric(readline(prompt = "Type I control-mean difference [1.0]: "))
  if (is.na(discordance)) discordance <- 1.0
  n_sim <- as.integer(readline(prompt = "MC iterations [2000]: "))
  if (is.na(n_sim)) n_sim <- 2000
  list(ratio = ratio, n_total = n_total, n_conventional = n_conventional,
       sd_int = sd_int, delta_tx = delta_tx, discordance = discordance,
       exclude = exclude, method = method,
       n_sim = n_sim, seed = 42, export = NULL, interactive = TRUE, compare = TRUE)
}

# ---- 6. main ----
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- parse_cli(args)
  if (opts$help) {
    cat(readLines(file.path(.sd, "borrowing_calculator.R"))[2:22], sep = "\n")
    return(invisible())
  }
  params <- if (opts$interactive) interactive_menu() else opts
  res <- do.call(calculate_borrowing, params[intersect(
    names(params),
    c("ratio", "n_total", "n_conventional", "sd_int", "delta_tx", "discordance",
      "exclude", "method", "n_sim", "seed", "compare"))])
  print_results(res)
  if (!is.null(params$export)) {
    row <- data.frame(
      ratio = res$ratio_label, n_total = res$n_total,
      n_conventional = res$n_conventional, n_placebo = res$n_placebo,
      sd_int = res$sd_int, delta_tx = res$delta_tx, discordance = res$discordance,
      pool_max_departure = res$pool_max_departure, method = res$method,
      meta_mu = res$meta_mu, meta_tau2 = res$meta_tau2, meta_I2 = res$meta_I2,
      power = res$power, no_borrow_power = res$no_borrow_power,
      trad_power = res$trad_power, borrowing_gain = res$borrowing_gain,
      placebo_reduction = res$placebo_reduction,
      type1_disc0 = res$type1_disc0, type1_disc1 = res$type1_disc1)
    write.csv(row, params$export, row.names = FALSE)
    cat("Saved main result row to:", params$export, "\n")
  }
}

if (sys.nframe() == 0) main()
