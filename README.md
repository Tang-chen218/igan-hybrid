# IgAN-Hybrid

Decision-support calculator for hybrid IgA nephropathy (IgAN) trial designs that
augment the concurrent control arm with external control data via elastic-prior
borrowing (Jiang, Nie, Yuan, *Biometrics* 2023).

Companion tool to the manuscript:

> *Reducing Placebo Enrollment in IgA Nephropathy Trials by Augmenting the Concurrent Control With External Control Data*

Given sample size, allocation ratio, treatment effect, within-trial SD, and an
external control pool, the tool returns power, type I error, and placebo reduction.
It is a planning aid, not a guarantee of any fixed power target.

## Requirements

- R (>= 4.0)
- CLI / functions: base R only
- Browser UI: [`shiny`](https://cran.r-project.org/package=shiny)

## Files

| File | Purpose |
|------|---------|
| `app.R` | Local Shiny UI |
| `borrowing_calculator.R` | CLI, interactive menu, and `calculate_borrowing()` |
| `pool.R` | Seven-study external control pool and random-effects meta-analysis |
| `simulate.R` | Monte Carlo elastic-prior engine |

## Quick start (browser)

```bash
cd igan-hybrid
# install.packages("shiny")   # once
Rscript -e 'shiny::runApp(".", launch.browser = TRUE)'
```

Or open `app.R` in RStudio and click **Run App**. Choose inputs, then **Run**.

## Quick start (CLI)

```bash
Rscript borrowing_calculator.R
Rscript borrowing_calculator.R --ratio 2to1 --n-total 300 --sd-int 6.0 \
  --delta-tx 2.0 --discordance 1.0 --n-sim 10000
Rscript borrowing_calculator.R --exclude APPLAUSE-IgAN --ratio 1to1 --n-sim 10000
Rscript borrowing_calculator.R --method none
Rscript borrowing_calculator.R --interactive
Rscript borrowing_calculator.R --help
```

### R function

```r
source("borrowing_calculator.R")
r <- calculate_borrowing(
  ratio = "1to1", n_total = 300, n_conventional = 370,
  sd_int = 6.0, delta_tx = 2.0, discordance = 1.0,
  n_sim = 10000, compare = FALSE
)
r$power
r$type1_disc0
r$type1_disc1
r$placebo_reduction
```

## Outputs

- External pool meta-analysis (μ, τ², I²)
- Hybrid power (Monte Carlo); no-borrowing and conventional 1:1 power (analytic)
- Placebo reduction versus the conventional design
- Type I error when control means are equal and at a user-set mean difference
- Power versus SD_internal (5.0–8.0) in the Shiny UI
- Optional 1:1 / 2:1 / 3:1 comparison at fixed total n

## Verification (n_sim = 10000, seed = 42)

| Config | Calculator | Manuscript |
|--------|-----------:|-----------:|
| 1:1, full completed-trials pool | 85.1% | 85.1% |
| 2:1, full completed-trials pool | 83.0% | 83.0% |
| 1:1, excluding APPLAUSE-IgAN | 89.6% / type I 5.9% | 89.6% / 5.9% |

## Status

- CLI and local Shiny UI: available
- Code: https://github.com/Tang-chen218/igan-hybrid
- Hosted app (shinyapps account `pkufh-igan`, app name **`igan-hybrid-calculator`**):

```r
rsconnect::deployApp(
  appDir = ".",
  appName = "igan-hybrid-calculator",
  appTitle = "IgAN-Hybrid Calculator",
  forceUpdate = TRUE
)
```

Expected URL: `https://pkufh-igan.shinyapps.io/igan-hybrid-calculator/`
