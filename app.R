# IgAN-Hybrid — 本地浏览器界面（Phase 1B）
# 启动：Rscript -e 'shiny::runApp(".", launch.browser=TRUE)'
library(shiny)
library(shinyglass)

source("borrowing_calculator.R", local = FALSE)

pct <- function(x) if (is.na(x)) "—" else sprintf("%.1f%%", 100 * x)
pp <- function(x) if (is.na(x)) "—" else sprintf("%+.1f percentage points", 100 * x)

# 组成与 Supplementary Table S5 脚注一致；去掉 APPLAUSE 用 Custom
POOL_PRESETS <- list(
  seven = EXT_STUDIES$author,
  chronic = c("TESTING", "NefIgArd", "PROTECT", "EMPA-KIDNEY IgAN"),
  total = c("DAPA-CKD IgAN", "APPLAUSE-IgAN", "ALIGN"),
  dedicated = setdiff(EXT_STUDIES$author, c("DAPA-CKD IgAN", "EMPA-KIDNEY IgAN"))
)

# 决策主指标：大数字；次要用紧凑对照表，避免 power 相关条目重复罗列
verdict_row <- function(label, value) {
  tags$div(
    style = "display:flex; justify-content:space-between; align-items:baseline; margin:8px 0;",
    tags$span(style = "font-size:1.05em;", label),
    tags$span(style = "font-size:1.35em; font-weight:600;", value)
  )
}

ui <- fluidPage(
  theme = glass_theme(preset = "light", intensity = 0.35),
  titlePanel("IgAN-Hybrid Calculator"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "pool_preset", "External control pool",
        choices = c(
          "All completed trials" = "seven",
          "Chronic eGFR-slope trials only" = "chronic",
          "Total eGFR-slope trials only" = "total",
          "Dedicated IgAN trials (excluding CKD subgroups)" = "dedicated",
          "Custom selection" = "custom"
        ),
        selected = "seven"
      ),
      checkboxGroupInput(
        "keep", "Trials included",
        choices = setNames(EXT_STUDIES$author, EXT_STUDIES$author),
        selected = EXT_STUDIES$author
      ),
      helpText("Custom: uncheck any trial to drop it."),
      radioButtons("ratio", "Randomization ratio (treatment:control)",
                   choices = c("1:1" = "1to1", "2:1" = "2to1", "3:1" = "3to1"),
                   selected = "1to1", inline = TRUE),
      numericInput("n_conventional", "Total sample size (conventional 1:1)",
                   value = 370, min = 50, max = 2000, step = 10),
      numericInput("n_total", "Total sample size (hybrid trial)", value = 300, min = 50, max = 1000, step = 10),
      numericInput("sd_int", "Planning SD_internal (for operating characteristics)",
                   value = 6.0, min = 3, max = 12, step = 0.5),
      helpText("Point estimate and type I error use this SD. The power curve below shows sensitivity across 5.0–8.0."),
      numericInput("delta_tx", "Assumed treatment effect on annualized eGFR slope",
                   value = 2.0, min = 0, max = 5, step = 0.5),
      numericInput("discordance", "Control-mean difference for type I error",
                   value = 1.0, min = 0, max = 5, step = 0.1),
      uiOutput("discordance_hint"),
      helpText("Effect and mean difference: mL/min/1.73 m²/yr."),
      radioButtons("method", "Borrowing method",
                   choices = c("Elastic prior (hybrid)" = "elastic",
                               "No borrowing (analytic power only)" = "none"),
                   selected = "elastic"),
      numericInput("n_sim", "Monte Carlo iterations", value = 10000, min = 200, max = 10000, step = 500),
      checkboxInput("compare", "Also compare 1:1, 2:1, and 3:1 at this sample size", value = FALSE),
      actionButton("go", "Run", class = "btn-primary"),
      width = 4
    ),
    mainPanel(
      width = 8,
      tags$p(tags$em("Estimates under the selected design assumptions.")),
      uiOutput("results_ui"),
      uiOutput("sd_curve_header"),
      plotOutput("sd_curve_plot", height = "340px"),
      uiOutput("cmp_header"),
      tableOutput("cmp_table")
    )
  )
)

# 全文池用稿中 1.0；其他池用该池 max|yi−μ|（一位小数）
suggested_discordance <- function(keep) {
  if (is.null(keep) || length(keep) < 2L) return(1.0)
  if (setequal(keep, EXT_STUDIES$author)) return(1.0)
  st <- EXT_STUDIES[EXT_STUDIES$author %in% keep, , drop = FALSE]
  meta <- meta_for_studies(st)
  round(pool_max_abs_departure(st, meta$mu_re), 1)
}

server <- function(input, output, session) {
  observeEvent(input$pool_preset, {
    if (input$pool_preset == "custom") return()
    updateCheckboxGroupInput(session, "keep", selected = POOL_PRESETS[[input$pool_preset]])
  }, ignoreInit = FALSE)

  observeEvent(input$keep, {
    if (is.null(input$keep)) return()
    preset <- input$pool_preset
    if (preset != "custom") {
      target <- POOL_PRESETS[[preset]]
      if (!setequal(input$keep, target)) {
        updateSelectInput(session, "pool_preset", selected = "custom")
      }
    }
    # 换池自动对齐 discordance，结果标题才会跟当前池走
    updateNumericInput(session, "discordance", value = suggested_discordance(input$keep))
  }, ignoreInit = FALSE)

  output$discordance_hint <- renderUI({
    keep <- input$keep
    if (is.null(keep) || length(keep) < 2) {
      return(helpText("Select ≥2 trials to set a pool-specific mean difference."))
    }
    st <- EXT_STUDIES[EXT_STUDIES$author %in% keep, , drop = FALSE]
    meta <- meta_for_studies(st)
    mx <- pool_max_abs_departure(st, meta$mu_re)
    sug <- suggested_discordance(keep)
    if (setequal(keep, EXT_STUDIES$author)) {
      helpText(sprintf(
        "Largest |trial slope − pooled mean| = %.2f; default 1.0 matches the manuscript.",
        mx
      ))
    } else {
      helpText(sprintf(
        "Largest |trial slope − pooled mean| = %.2f; default set to %.1f for this pool.",
        mx, sug
      ))
    }
  })

  res <- eventReactive(input$go, {
    keep <- input$keep
    if (is.null(keep) || length(keep) < 2) {
      validate(need(FALSE, "Keep at least 2 studies in the pool."))
    }
    exclude <- setdiff(EXT_STUDIES$author, keep)
    withProgress(message = "Running simulations…", value = 0.2, {
      out <- calculate_borrowing(
        ratio = input$ratio,
        n_total = as.integer(input$n_total),
        n_conventional = as.integer(input$n_conventional),
        sd_int = as.numeric(input$sd_int),
        delta_tx = as.numeric(input$delta_tx),
        discordance = as.numeric(input$discordance),
        exclude = exclude,
        method = input$method,
        n_sim = as.integer(input$n_sim),
        seed = 42,
        compare = isTRUE(input$compare) && input$method == "elastic"
      )
      setProgress(1)
      out
    })
  })

  output$results_ui <- renderUI({
    r <- res()
    method_lab <- if (r$method == "elastic") {
      "Elastic prior borrowing from the external control pool."
    } else {
      "No external borrowing; analytic power only."
    }

    # 一条对照句代替两个独立 delta 指标
    cmp_note <- if (r$method == "none" || is.na(r$borrowing_gain)) {
      NULL
    } else {
      sprintf(
        "Compared with the same total sample size without external borrowing: %s. Compared with the conventional 1:1 design (n = %d): %s.",
        pp(r$borrowing_gain), r$n_conventional, pp(r$power_vs_trad)
      )
    }

    type1_pair <- if (r$method == "none") {
      "—"
    } else {
      sprintf("%s / %s", pct(r$type1_disc0), pct(r$type1_disc1))
    }

    tagList(
      tags$h3("1. External control pool"),
      tags$p(sprintf(
        "%d trial(s); %s participants in external control arms. Random-effects pooled eGFR slope %.2f mL/min/1.73 m²/yr (τ² = %.3f; I² = %.1f%%).",
        r$pool_n_studies, format(r$pool_n, big.mark = ","),
        r$meta_mu, r$meta_tau2, r$meta_I2
      )),
      tags$p(sprintf(
        "Hybrid design: %s, n = %d (%d treatment : %d concurrent control). Conventional reference: 1:1, n = %d. Assumed effect = %.1f; SD_internal = %.1f.",
        r$ratio_label, r$n_total, r$n_treatment, r$n_placebo,
        r$n_conventional, r$delta_tx, r$sd_int
      )),
      tags$p(method_lab),

      tags$h3("2. Operating characteristics"),
      verdict_row("Power with external borrowing", pct(r$power)),
      verdict_row(
        sprintf("Reduction in concurrent placebo vs conventional n = %d", r$n_conventional),
        sprintf("%d%% fewer (%d vs %d)",
                round(100 * r$placebo_reduction), r$n_placebo, r$n_conventional_placebo)
      ),
      verdict_row(
        sprintf("Type I error (control means equal / differ by %.1f)", r$discordance),
        type1_pair
      ),

      tags$h3("3. Power relative to reference designs"),
      tags$table(
        style = "width:100%; border-collapse:collapse; margin-bottom:6px;",
        tags$thead(tags$tr(
          tags$th(style = "text-align:left; padding:6px 0; border-bottom:1px solid #ddd;", "Design"),
          tags$th(style = "text-align:right; padding:6px 0; border-bottom:1px solid #ddd;", "Power")
        )),
        tags$tbody(
          tags$tr(
            tags$td(style = "padding:6px 0;",
                    sprintf("Hybrid %s, n = %d (with borrowing)", r$ratio_label, r$n_total)),
            tags$td(style = "text-align:right; padding:6px 0;", pct(r$power))
          ),
          tags$tr(
            tags$td(style = "padding:6px 0; color:#555;",
                    sprintf("Same n = %d, without external borrowing", r$n_total)),
            tags$td(style = "text-align:right; padding:6px 0; color:#555;", pct(r$no_borrow_power))
          ),
          tags$tr(
            tags$td(style = "padding:6px 0; color:#555;",
                    sprintf("Conventional 1:1, n = %d", r$n_conventional)),
            tags$td(style = "text-align:right; padding:6px 0; color:#555;", pct(r$trad_power))
          )
        )
      ),
      if (!is.null(cmp_note)) {
        tags$p(style = "color:#555; font-size:0.9em;", cmp_note)
      }
    )
  })

  output$sd_curve_header <- renderUI({
    r <- res()
    req(r)
    if (is.null(r$sd_curve)) return(NULL)
    tagList(
      tags$h3("4. Power versus SD_internal"),
      tags$p(
        style = "color:#555; font-size:0.9em;",
        "Power as a function of SD_internal (5.0–8.0)."
      )
    )
  })

  output$sd_curve_plot <- renderPlot({
    r <- res()
    req(r)
    d <- r$sd_curve
    req(d)

    # 与 Figure S2 一致：黑=conventional，橙=无借用，蓝=借用
    col_trad <- "#000000"
    col_nob  <- "#E69F00"
    col_hyb  <- "#0072B2"
    sd_sel <- r$sd_int

    op <- par(mar = c(4.2, 4.2, 1.2, 1.2), family = "sans")
    on.exit(par(op), add = TRUE)
    plot(d$sd_int, d$conventional, type = "n",
         xlim = c(4.8, 8.2), ylim = c(0.40, 1.00),
         xlab = expression(SD[internal]~"(mL/min/1.73 m"^2*"/yr)"),
         ylab = "Power", las = 1)
    abline(h = 0.80, lty = 2, col = "grey45")
    abline(v = sd_sel, lty = 3, col = "grey60")
    lines(d$sd_int, d$conventional, col = col_trad, lwd = 2)
    points(d$sd_int, d$conventional, pch = 21, bg = "white", col = col_trad, cex = 1.1)
    lines(d$sd_int, d$no_borrow, col = col_nob, lwd = 2, lty = 2)
    points(d$sd_int, d$no_borrow, pch = 21, bg = "white", col = col_nob, cex = 1.1)
    lines(d$sd_int, d$hybrid, col = col_hyb, lwd = 2)
    points(d$sd_int, d$hybrid, pch = 21, bg = "white", col = col_hyb, cex = 1.1)
    text(7.95, 0.835, "80%", col = "grey45", cex = 0.85, adj = 1)
    legend("bottomleft",
           legend = c(
             sprintf("Conventional 1:1, n = %d", r$n_conventional),
             sprintf("Without borrowing, n = %d", r$n_total),
             sprintf("With borrowing, n = %d", r$n_total)
           ),
           col = c(col_trad, col_nob, col_hyb), lwd = 2, lty = c(1, 2, 1),
           pch = 21, pt.bg = "white", bty = "n", cex = 0.85)
  })

  output$cmp_header <- renderUI({
    r <- res()
    if (is.null(r$comparison)) return(NULL)
    tagList(
      tags$h3("5. Allocation ratio comparison"),
      tags$p(style = "color:#555; font-size:0.92em;",
             "Power with borrowing at fixed total n across randomization ratios.")
    )
  })

  output$cmp_table <- renderTable({
    r <- res()
    if (is.null(r$comparison)) return(NULL)
    data.frame(
      `Randomization ratio` = r$comparison$ratio,
      `Power with borrowing` = sprintf("%.1f%%", 100 * r$comparison$power),
      `Concurrent placebo (n)` = r$comparison$placebo_n,
      `Power vs conventional` = sprintf("%+.1f pp", 100 * r$comparison$power_vs_trad),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, width = "100%")
}

shinyApp(ui, server)
