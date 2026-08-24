#!/usr/bin/env Rscript
# 部署到 shinyapps.io（账户 pkufh-igan，应用 igan-hybrid-calculator）
# 1) 在 shinyapps.io → Tokens → Show，复制 setAccountInfo(...) 到 R 里执行一次
# 2) 再运行：Rscript deploy_shinyapps.R

stopifnot(requireNamespace("rsconnect", quietly = TRUE))
stopifnot(requireNamespace("shiny", quietly = TRUE))

accts <- rsconnect::accounts()
if (nrow(accts) == 0) {
  stop(
    "No shinyapps account configured.\n",
    "Run the setAccountInfo(...) line from shinyapps.io → Tokens in R once, then retry."
  )
}

app_dir <- normalizePath(".", winslash = "/")
rsconnect::deployApp(
  appDir = app_dir,
  appName = "igan-hybrid-calculator",
  appTitle = "IgAN-Hybrid Calculator",
  forceUpdate = TRUE
)
