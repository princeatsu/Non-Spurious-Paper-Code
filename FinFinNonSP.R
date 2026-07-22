################################################################################
#  FULL ANALYSIS SCRIPT
#  Tables 1-4  : Item-level descriptives + all estimators (single dataset)
#  Tables 5-6  : Mean Bias and RMSE from simulation
#  Figures 1-2 : Bias and RMSE line charts
#
#  KEY CHANGES FROM PREVIOUS VERSION:
#  (1) IC CALIBRATION: loadings are iteratively adjusted before any analysis
#      so the achieved Cronbach alpha is approximately 0.80 for every test
#      length, matching the claim in Section 2.4 of the manuscript.
#  (2) CONSISTENT HEADERS: all six tables use the same column labels:
#      r-raw | Guilford | Cureton | Henrysson | CITC | AIOSR
#  (3) APA FORMATTING: Tables 5 and 6 now use identical APA styling
#      (borders, fonts, column widths) as Tables 1-4.
################################################################################


################################
# 1. PACKAGES
################################

required <- c("psych", "flextable", "officer", "ggplot2", "tidyr", "dplyr", "scales")
new_pkgs  <- required[!required %in% rownames(installed.packages())]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
suppressPackageStartupMessages(
  invisible(lapply(required, library, character.only = TRUE))
)


################################
# 2. USER SETTINGS
################################

N            <- 1000L
n_items_grid <- c(5L, 10L, 20L, 30L)
seed_base    <- 2026L
R_sim        <- 500L          # simulation replications (use 100 for quick test)
target_alpha <- 0.80          # desired internal consistency

out_dir       <- "C:/Users/Akpabli Prince Atsu/Desktop/Prof Gabriel Asare Okyere/Articles/Non_Spurious/Psychometrika"
out_file_t14  <- file.path(out_dir, "Tables_1_to_4_APA.docx")
out_file_t56  <- file.path(out_dir, "Tables_Bias_RMSE_APA.docx")
fig_bias_path <- file.path(out_dir, "Figure_Bias.png")
fig_rmse_path <- file.path(out_dir, "Figure_RMSE.png")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


################################
# 3. CORRECTION FORMULAE
################################

# Equation 1 - Guilford (1953)
guilford_fn <- function(r_it_v, sigma_i_v, sigma_t) {
  num <- sigma_t * r_it_v - sigma_i_v
  den <- sqrt(sigma_t^2 + sigma_i_v^2 - 2 * sigma_t * sigma_i_v * r_it_v)
  ifelse(den < 1e-12, NA_real_, num / den)
}

# Equation 2 - Henrysson (1963)
henrysson_fn <- function(r_it_v, sigma_i_v, sigma_t) {
  n_i <- length(sigma_i_v)
  num <- sigma_t * r_it_v - sigma_i_v
  den <- sqrt(sigma_t^2 - sum(sigma_i_v^2))
  if (!is.finite(den) || den < 1e-12) return(rep(NA_real_, n_i))
  num / den * sqrt(n_i / (n_i - 1L))
}

# Equation 3 - Cureton (1966)
cureton_fn <- function(r_it_v, sigma_i_v, sigma_t, r_tt) {
  lam  <- sigma_i_v / sigma_t
  disc <- r_tt^2 - 4 * r_tt * lam * (r_it_v - lam)
  disc <- pmax(disc, 0)
  ifelse(lam < 1e-12, NA_real_, (r_tt - sqrt(disc)) / (2 * lam))
}

# Equation 9 - internal-consistency formula
rtt_eq9 <- function(r_it_v, sigma_i_v, sigma_t) {
  denom <- sigma_t^2 - sum(r_it_v^2 * sigma_i_v^2)
  if (!is.finite(denom) || abs(denom) < 1e-12) return(NA_real_)
  (sigma_t^2 - sum(sigma_i_v^2)) / denom
}

# Equation 11 - CITC
citc_fn <- function(r_it_v, sigma_i_v, sigma_t) {
  r_tt   <- rtt_eq9(r_it_v, sigma_i_v, sigma_t)
  if (!is.finite(r_tt)) return(rep(NA_real_, length(r_it_v)))
  r_ii_v <- r_it_v^2 * r_tt
  r_it_v - (sigma_i_v / sigma_t) * (1 - r_ii_v)
}

# Equation 12 - AIOSR
aiosr_fn <- function(r_it_v, sigma_i_v, sigma_t) {
  r_tt <- rtt_eq9(r_it_v, sigma_i_v, sigma_t)
  if (!is.finite(r_tt) || r_tt <= 0) return(rep(NA_real_, length(r_it_v)))
  citc_fn(r_it_v, sigma_i_v, sigma_t) / sqrt(r_tt)
}


################################
# 4. IC CALIBRATION
#
#    Runs n_calib quick replications (N_calib examinees each) and
#    adjusts the base factor loading until the mean Cronbach alpha
#    is within tol of target_alpha. This ensures that Section 2.4's
#    claim "internal consistency approximately 0.80" is enforced by
#    the code, not just intended.
################################

target_loading_analytical <- function(k) sqrt(target_alpha / (target_alpha * (k - 1) / k + 1 / k))

calibrate_loading <- function(n_items,
                              N_calib  = 500L,
                              n_calib  = 50L,
                              tol      = 0.015,
                              max_iter = 15L,
                              calib_seed = 42L) {
  lam <- target_loading_analytical(n_items)
  message(sprintf("\nCalibrating loading for n = %d items (target alpha = %.2f)...",
                  n_items, target_alpha))
  
  mean_alpha <- NA_real_
  
  for (iter in seq_len(max_iter)) {
    lam_iter <- lam
    
    alphas <- vapply(seq_len(n_calib), function(r) {
      set.seed(calib_seed + r + n_items * 100L)
      lam_j   <- pmin(pmax(stats::rnorm(n_items, lam_iter, 0.04), 0.20), 0.90)
      mu_j    <- stats::runif(n_items, 2, 5)
      sd_j    <- stats::runif(n_items, 1.25, 1.75)
      F_lat   <- stats::rnorm(N_calib)
      X <- vapply(seq_len(n_items), function(j) {
        eps <- stats::rnorm(N_calib, 0, sqrt(max(1 - lam_j[j]^2, 1e-6)))
        Z   <- mu_j[j] + sd_j[j] * (lam_j[j] * F_lat + eps)
        pmin(pmax(round(Z), 0L), 6L)
      }, numeric(N_calib))
      tryCatch(
        psych::alpha(as.data.frame(X),
                     check.keys = FALSE,
                     warnings   = FALSE)$total$raw_alpha,
        error = function(e) NA_real_
      )
    }, numeric(1))
    
    mean_alpha <- mean(alphas, na.rm = TRUE)
    message(sprintf("  Iter %2d: base_lam = %.4f,  mean alpha = %.3f",
                    iter, lam_iter, mean_alpha))
    
    if (abs(mean_alpha - target_alpha) <= tol) {
      message(sprintf("  Converged. Final base_lam = %.4f, expected alpha = %.3f\n",
                      lam, mean_alpha))
      return(lam)
    }
    
    # Proportional correction: scale loading up/down to move alpha toward target
    lam <- pmin(pmax(lam * sqrt(target_alpha / mean_alpha), 0.20), 0.90)
  }
  
  warning(sprintf(
    "n = %d: calibration did not converge after %d iterations. Using lam = %.4f (mean alpha = %.3f).",
    n_items, max_iter, lam, mean_alpha))
  return(lam)
}

# Run calibration once for each test length and store results
message("=== PHASE 1: IC CALIBRATION ===")
calibrated_lams <- setNames(
  vapply(n_items_grid, calibrate_loading, numeric(1)),
  as.character(n_items_grid)
)
message("Calibration complete. Calibrated base loadings:")
print(round(calibrated_lams, 4))


################################
# 5. DATA GENERATION
#    Accepts pre-calibrated base loading so all datasets share
#    the same target IC.
################################

generate_data <- function(N, n_items, seed = NULL,
                          base_lam = calibrated_lams[as.character(n_items)]) {
  if (!is.null(seed)) set.seed(seed)
  
  loadings <- pmin(pmax(stats::rnorm(n_items, base_lam, 0.04), 0.20), 0.90)
  mu_j     <- stats::runif(n_items, min = 2,    max = 5)
  sd_j     <- stats::runif(n_items, min = 1.25, max = 1.75)
  F_latent <- stats::rnorm(N)
  
  X <- vapply(seq_len(n_items), function(j) {
    eps <- stats::rnorm(N, mean = 0,
                        sd   = sqrt(max(1 - loadings[j]^2, 1e-6)))
    Z   <- mu_j[j] + sd_j[j] * (loadings[j] * F_latent + eps)
    pmin(pmax(round(Z), 0L), 6L)
  }, numeric(N))
  
  colnames(X) <- paste0("q", seq_len(n_items))
  list(df = as.data.frame(X), F_latent = F_latent, loadings = loadings)
}


################################
# 6. COMPUTE ALL ESTIMATORS
################################

compute_estimators <- function(df, alpha_tt) {
  Y        <- rowSums(df)
  sigma_t  <- stats::sd(Y)
  sigma_iv <- vapply(df, stats::sd, numeric(1))
  r_it_v   <- vapply(df, function(x) stats::cor(x, Y), numeric(1))
  
  data.frame(
    item      = names(df),
    r_it      = r_it_v,
    Guilford  = guilford_fn(r_it_v, sigma_iv, sigma_t),
    Cureton   = cureton_fn(r_it_v, sigma_iv, sigma_t, alpha_tt),
    Henrysson = henrysson_fn(r_it_v, sigma_iv, sigma_t),
    CITC      = citc_fn(r_it_v, sigma_iv, sigma_t),
    AIOSR     = aiosr_fn(r_it_v, sigma_iv, sigma_t),
    M         = colMeans(df),
    SD        = sigma_iv,
    n         = nrow(df),
    stringsAsFactors = FALSE
  )
}


################################
# 7. TABLES 1-4
#    Single fixed dataset per test length.
#    Column headers now match the Bias/RMSE tables:
#    r-raw | Guilford | Cureton | Henrysson | CITC | AIOSR
#
#    NOTE FOR MANUSCRIPT: the column key in Section 3 Results
#    should be updated to reflect these new header names.
################################

message("\n=== PHASE 2: GENERATING TABLES 1-4 ===")

results_list <- lapply(n_items_grid, function(k) {
  gen      <- generate_data(N = N, n_items = k, seed = seed_base + k)
  alpha_tt <- tryCatch(
    psych::alpha(gen$df, check.keys = FALSE, warnings = FALSE)$total$raw_alpha,
    error = function(e) NA_real_
  )
  est <- compute_estimators(gen$df, alpha_tt)
  
  tbl <- data.frame(
    "Item"      = est$item,
    "n"         = est$n,
    "M"         = round(est$M,        2),
    "SD"        = round(est$SD,       2),
    "r-raw"     = round(est$r_it,     3),
    "Guilford"  = round(est$Guilford, 3),
    "Cureton"   = round(est$Cureton,  3),
    "Henrysson" = round(est$Henrysson,3),
    "CITC"      = round(est$CITC,     3),
    "AIOSR"     = round(est$AIOSR,    3),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  list(tbl = tbl, alpha = round(alpha_tt, 2), k = k)
})
names(results_list) <- as.character(n_items_grid)


################################
# 8. SIMULATION LOOP
#    Same generate_data() and calibrated loadings as Tables 1-4,
#    so simulation results are fully coherent with the main tables.
################################

message(sprintf("\n=== PHASE 3: SIMULATION (%d replications per condition) ===", R_sim))

sim_summary <- lapply(n_items_grid, function(k) {
  message(sprintf("  n_items = %d", k))
  
  rep_rows <- vector("list", R_sim)
  
  for (r in seq_len(R_sim)) {
    gen      <- generate_data(N = N, n_items = k,
                              seed = seed_base + k * 1000L + r)
    alpha_tt <- tryCatch(
      psych::alpha(gen$df, check.keys = FALSE, warnings = FALSE)$total$raw_alpha,
      error = function(e) NA_real_
    )
    est       <- compute_estimators(gen$df, alpha_tt)
    true_corr <- vapply(gen$df,
                        function(x) stats::cor(x, gen$F_latent),
                        numeric(1))
    
    rep_rows[[r]] <- data.frame(
      rep       = r,
      item      = est$item,
      true_corr = true_corr,
      r_it      = est$r_it,
      Guilford  = est$Guilford,
      Cureton   = est$Cureton,
      Henrysson = est$Henrysson,
      CITC      = est$CITC,
      AIOSR     = est$AIOSR,
      stringsAsFactors = FALSE
    )
  }
  
  df_all <- do.call(rbind, rep_rows)
  
  per_item <- lapply(split(df_all, df_all$item), function(d) {
    tc <- d$true_corr
    data.frame(
      item           = d$item[1],
      bias_r_it      = mean(d$r_it      - tc, na.rm = TRUE),
      bias_Guilford  = mean(d$Guilford  - tc, na.rm = TRUE),
      bias_Cureton   = mean(d$Cureton   - tc, na.rm = TRUE),
      bias_Henrysson = mean(d$Henrysson - tc, na.rm = TRUE),
      bias_CITC      = mean(d$CITC      - tc, na.rm = TRUE),
      bias_AIOSR     = mean(d$AIOSR     - tc, na.rm = TRUE),
      rmse_r_it      = sqrt(mean((d$r_it      - tc)^2, na.rm = TRUE)),
      rmse_Guilford  = sqrt(mean((d$Guilford  - tc)^2, na.rm = TRUE)),
      rmse_Cureton   = sqrt(mean((d$Cureton   - tc)^2, na.rm = TRUE)),
      rmse_Henrysson = sqrt(mean((d$Henrysson - tc)^2, na.rm = TRUE)),
      rmse_CITC      = sqrt(mean((d$CITC      - tc)^2, na.rm = TRUE)),
      rmse_AIOSR     = sqrt(mean((d$AIOSR     - tc)^2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })
  per_item_df <- do.call(rbind, per_item)
  
  data.frame(
    n_items        = k,
    bias_r_it      = round(mean(per_item_df$bias_r_it),      3),
    bias_Guilford  = round(mean(per_item_df$bias_Guilford),   3),
    bias_Cureton   = round(mean(per_item_df$bias_Cureton),    3),
    bias_Henrysson = round(mean(per_item_df$bias_Henrysson),  3),
    bias_CITC      = round(mean(per_item_df$bias_CITC),       3),
    bias_AIOSR     = round(mean(per_item_df$bias_AIOSR),      3),
    rmse_r_it      = round(mean(per_item_df$rmse_r_it),       3),
    rmse_Guilford  = round(mean(per_item_df$rmse_Guilford),   3),
    rmse_Cureton   = round(mean(per_item_df$rmse_Cureton),    3),
    rmse_Henrysson = round(mean(per_item_df$rmse_Henrysson),  3),
    rmse_CITC      = round(mean(per_item_df$rmse_CITC),       3),
    rmse_AIOSR     = round(mean(per_item_df$rmse_AIOSR),      3),
    stringsAsFactors = FALSE
  )
})

aggregate_summaries <- do.call(rbind, sim_summary)
message("Simulation complete.")

write.csv(aggregate_summaries,
          file.path(out_dir, "aggregate_summary_by_nitems.csv"),
          row.names = FALSE)


################################
# 9. BIAS AND RMSE DISPLAY TABLES
#    Column order matches Tables 1-4: Guilford | Cureton | Henrysson | CITC | AIOSR
################################

bias_tbl <- data.frame(
  "k"         = aggregate_summaries$n_items,
  "r-raw"     = aggregate_summaries$bias_r_it,
  "Guilford"  = aggregate_summaries$bias_Guilford,
  "Cureton"   = aggregate_summaries$bias_Cureton,
  "Henrysson" = aggregate_summaries$bias_Henrysson,
  "CITC"      = aggregate_summaries$bias_CITC,
  "AIOSR"     = aggregate_summaries$bias_AIOSR,
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

rmse_tbl <- data.frame(
  "k"         = aggregate_summaries$n_items,
  "r-raw"     = aggregate_summaries$rmse_r_it,
  "Guilford"  = aggregate_summaries$rmse_Guilford,
  "Cureton"   = aggregate_summaries$rmse_Cureton,
  "Henrysson" = aggregate_summaries$rmse_Henrysson,
  "CITC"      = aggregate_summaries$rmse_CITC,
  "AIOSR"     = aggregate_summaries$rmse_AIOSR,
  check.names      = FALSE,
  stringsAsFactors = FALSE
)


################################
# 10. FIGURES
#     Line plots with colour-blind-friendly palette,
#     distinct line types and point shapes for greyscale printing.
################################

estimator_order  <- c("r_it", "Guilford", "Cureton", "Henrysson", "CITC", "AIOSR")
estimator_labels <- c("r-raw", "Guilford", "Cureton", "Henrysson", "CITC", "AIOSR")

cbf_palette <- c(
  "r_it"      = "#000000",
  "Guilford"  = "#E69F00",
  "Cureton"   = "#009E73",
  "Henrysson" = "#56B4E9",
  "CITC"      = "#CC79A7",
  "AIOSR"     = "#FF0000"
)

line_types <- c(
  "r_it"      = "solid",
  "Guilford"  = "dashed",
  "Cureton"   = "longdash",
  "Henrysson" = "dotdash",
  "CITC"      = "twodash",
  "AIOSR"     = "solid"
)

point_shapes <- c(
  "r_it"      = 16,
  "Guilford"  = 17,
  "Cureton"   = 18,
  "Henrysson" = 15,
  "CITC"      = 8,
  "AIOSR"     = 21
)

apa_theme <- ggplot2::theme_classic(base_family = "serif", base_size = 12) +
  ggplot2::theme(
    axis.title         = ggplot2::element_text(size = 12),
    axis.text          = ggplot2::element_text(size = 11, color = "black"),
    legend.title       = ggplot2::element_blank(),
    legend.text        = ggplot2::element_text(size = 10),
    legend.position    = "bottom",
    legend.key.width   = ggplot2::unit(1.2, "cm"),
    legend.key.size    = ggplot2::unit(0.5, "cm"),
    panel.grid.major.y = ggplot2::element_line(color = "grey85", linewidth = 0.4),
    panel.grid.major.x = ggplot2::element_line(color = "grey92", linewidth = 0.3),
    plot.title         = ggplot2::element_blank()
  )

# Figure 1: Bias
bias_long <- tidyr::pivot_longer(
  aggregate_summaries,
  cols      = dplyr::starts_with("bias_"),
  names_to  = "estimator",
  values_to = "bias"
) %>%
  dplyr::mutate(
    estimator = sub("bias_", "", estimator),
    estimator = dplyr::recode(estimator, "r_it" = "r_it"),
    estimator = factor(estimator, levels = estimator_order)
  )

fig_bias <- ggplot2::ggplot(
  bias_long,
  ggplot2::aes(x = n_items, y = bias,
               colour   = estimator,
               linetype = estimator,
               shape    = estimator)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.5,
                      linetype = "solid", colour = "grey40") +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::geom_point(size = 2.8, fill = "white") +
  ggplot2::scale_colour_manual(values = cbf_palette, labels = estimator_labels) +
  ggplot2::scale_linetype_manual(values = line_types, labels = estimator_labels) +
  ggplot2::scale_shape_manual(values = point_shapes, labels = estimator_labels) +
  ggplot2::scale_x_continuous(breaks = n_items_grid,
                              labels = as.character(n_items_grid)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  ggplot2::labs(x = "Number of Items (k)",
                y = "Mean Bias (Estimator \u2212 True Correlation)") +
  apa_theme +
  ggplot2::guides(
    colour   = ggplot2::guide_legend(nrow = 2),
    linetype = ggplot2::guide_legend(nrow = 2),
    shape    = ggplot2::guide_legend(nrow = 2)
  )

ggplot2::ggsave(fig_bias_path, plot = fig_bias,
                width = 6.5, height = 4.5, dpi = 300, units = "in")
message("Saved: ", fig_bias_path)

# Figure 2: RMSE
rmse_long <- tidyr::pivot_longer(
  aggregate_summaries,
  cols      = dplyr::starts_with("rmse_"),
  names_to  = "estimator",
  values_to = "rmse"
) %>%
  dplyr::mutate(
    estimator = sub("rmse_", "", estimator),
    estimator = factor(estimator, levels = estimator_order)
  )

fig_rmse <- ggplot2::ggplot(
  rmse_long,
  ggplot2::aes(x = n_items, y = rmse,
               colour   = estimator,
               linetype = estimator,
               shape    = estimator)
) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::geom_point(size = 2.8, fill = "white") +
  ggplot2::scale_colour_manual(values = cbf_palette, labels = estimator_labels) +
  ggplot2::scale_linetype_manual(values = line_types, labels = estimator_labels) +
  ggplot2::scale_shape_manual(values = point_shapes, labels = estimator_labels) +
  ggplot2::scale_x_continuous(breaks = n_items_grid,
                              labels = as.character(n_items_grid)) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  ggplot2::labs(x = "Number of Items (k)",
                y = "Mean RMSE") +
  apa_theme +
  ggplot2::guides(
    colour   = ggplot2::guide_legend(nrow = 2),
    linetype = ggplot2::guide_legend(nrow = 2),
    shape    = ggplot2::guide_legend(nrow = 2)
  )

ggplot2::ggsave(fig_rmse_path, plot = fig_rmse,
                width = 6.5, height = 4.5, dpi = 300, units = "in")
message("Saved: ", fig_rmse_path)


################################
# 11. SHARED WORD HELPERS
################################

tnr <- "Times New Roman"

# APA flextable — Tables 1-4 (10 columns)
make_apa_ft_t14 <- function(df) {
  ft <- flextable::flextable(df)
  ft <- flextable::font(ft,       fontname = tnr, part = "all")
  ft <- flextable::fontsize(ft,   size = 11,      part = "all")
  ft <- flextable::italic(ft,     part = "header")
  ft <- flextable::bold(ft,       part = "header", bold = FALSE)
  ft <- flextable::align(ft,      align = "center", part = "all")
  ft <- flextable::align(ft, j = 1, align = "left", part = "all")
  ft <- flextable::padding(ft,    padding.top = 2, padding.bottom = 2, part = "all")
  ft <- flextable::border_remove(ft)
  thick <- officer::fp_border(color = "black", width = 1.5)
  thin  <- officer::fp_border(color = "black", width = 0.75)
  ft <- flextable::hline_top(ft,    part = "header", border = thick)
  ft <- flextable::hline_bottom(ft, part = "header", border = thin)
  ft <- flextable::hline_bottom(ft, part = "body",   border = thick)
  ft <- flextable::width(ft, j =  1, width = 0.55)   # Item
  ft <- flextable::width(ft, j =  2, width = 0.42)   # n
  ft <- flextable::width(ft, j =  3, width = 0.48)   # M
  ft <- flextable::width(ft, j =  4, width = 0.48)   # SD
  ft <- flextable::width(ft, j =  5, width = 0.62)   # r-raw
  ft <- flextable::width(ft, j =  6, width = 0.75)   # Guilford
  ft <- flextable::width(ft, j =  7, width = 0.70)   # Cureton
  ft <- flextable::width(ft, j =  8, width = 0.80)   # Henrysson
  ft <- flextable::width(ft, j =  9, width = 0.62)   # CITC
  ft <- flextable::width(ft, j = 10, width = 0.65)   # AIOSR
  ft
}

# APA flextable — Tables 5-6 (7 columns)
# Identical APA style to Tables 1-4 (same borders, fonts, sizes)
make_apa_ft_t56 <- function(df) {
  ft <- flextable::flextable(df)
  ft <- flextable::font(ft,       fontname = tnr, part = "all")
  ft <- flextable::fontsize(ft,   size = 11,      part = "all")
  ft <- flextable::italic(ft,     part = "header")
  ft <- flextable::bold(ft,       part = "header", bold = FALSE)
  ft <- flextable::align(ft,      align = "center", part = "all")
  ft <- flextable::align(ft, j = 1, align = "left", part = "all")
  ft <- flextable::padding(ft,    padding.top = 2, padding.bottom = 2, part = "all")
  ft <- flextable::border_remove(ft)
  thick <- officer::fp_border(color = "black", width = 1.5)
  thin  <- officer::fp_border(color = "black", width = 0.75)
  ft <- flextable::hline_top(ft,    part = "header", border = thick)
  ft <- flextable::hline_bottom(ft, part = "header", border = thin)
  ft <- flextable::hline_bottom(ft, part = "body",   border = thick)
  ft <- flextable::width(ft, j = 1, width = 0.55)   # k
  ft <- flextable::width(ft, j = 2, width = 0.80)   # r-raw
  ft <- flextable::width(ft, j = 3, width = 0.85)   # Guilford
  ft <- flextable::width(ft, j = 4, width = 0.82)   # Cureton
  ft <- flextable::width(ft, j = 5, width = 0.88)   # Henrysson
  ft <- flextable::width(ft, j = 6, width = 0.75)   # CITC
  ft <- flextable::width(ft, j = 7, width = 0.75)   # AIOSR
  ft
}

# Add one formatted paragraph
add_par <- function(doc, text, bold = FALSE, italic = FALSE,
                    size = 12, pad_after = 0, pad_before = 0) {
  officer::body_add_fpar(
    doc,
    officer::fpar(
      officer::ftext(text,
                     prop = officer::fp_text(bold        = bold,
                                             italic      = italic,
                                             font.size   = size,
                                             font.family = tnr)),
      fp_p = officer::fp_par(padding.bottom = pad_after,
                             padding.top    = pad_before,
                             line_spacing   = 1)
    )
  )
}

# Add APA Note paragraph (italic "Note. " + plain text)
add_note <- function(doc, note_text) {
  officer::body_add_fpar(
    doc,
    officer::fpar(
      officer::ftext("Note. ",
                     prop = officer::fp_text(italic      = TRUE,
                                             font.size   = 11,
                                             font.family = tnr)),
      officer::ftext(note_text,
                     prop = officer::fp_text(italic      = FALSE,
                                             font.size   = 11,
                                             font.family = tnr)),
      fp_p = officer::fp_par(padding.top = 4, line_spacing = 1)
    )
  )
}


################################
# 12. WORD DOCUMENT 1 — TABLES 1-4
################################

message("\n=== PHASE 4: BUILDING WORD DOCUMENTS ===")
message("Building Tables 1-4...")

doc1 <- officer::read_docx()

for (i in seq_along(n_items_grid)) {
  k   <- n_items_grid[i]
  obj <- results_list[[as.character(k)]]
  
  if (i > 1) doc1 <- officer::body_add_break(doc1)
  
  doc1 <- add_par(doc1, sprintf("Table %d", i), bold = TRUE, size = 12, pad_after = 0)
  doc1 <- add_par(doc1,
                  sprintf("Non-Corrected and Corrected Item\u2013Total Correlations for a Test of Length %d", k),
                  italic = TRUE, size = 12, pad_after = 6)
  doc1 <- flextable::body_add_flextable(doc1, value = make_apa_ft_t14(obj$tbl))
  doc1 <- add_note(doc1,
                   sprintf("Internal consistency of the test: Cronbach\u2019s \u03b1 = %.2f.", obj$alpha))
}

print(doc1, target = out_file_t14)
message(sprintf("Saved: %s", normalizePath(out_file_t14)))


################################
# 13. WORD DOCUMENT 2 — TABLES 5-6 + FIGURES 1-2
################################

message("Building Tables 5-6 and Figures...")

doc2 <- officer::read_docx()

# Table 5: Bias
doc2 <- add_par(doc2, "Table 5", bold = TRUE, size = 12, pad_after = 0)
doc2 <- add_par(doc2,
                paste0("Mean Bias of Item\u2013Total Correlation Estimators ",
                       "Across Test Lengths (Estimator \u2212 True Correlation)"),
                italic = TRUE, size = 12, pad_after = 6)
doc2 <- flextable::body_add_flextable(doc2, value = make_apa_ft_t56(bias_tbl))
doc2 <- add_note(doc2,
                 paste0("k = number of items. Values are mean bias averaged over items and ",
                        R_sim, " simulation replications (N = ", N, " per replication). ",
                        "Bias = estimator minus true item\u2013factor correlation; ",
                        "positive values indicate overestimation, negative values indicate underestimation. ",
                        "CITC = corrected item\u2013total correlation (Eq.\u00a011); ",
                        "AIOSR = adjusted for item overlap and scale reliability (Eq.\u00a012)."))

doc2 <- officer::body_add_break(doc2)

# Table 6: RMSE
doc2 <- add_par(doc2, "Table 6", bold = TRUE, size = 12, pad_after = 0)
doc2 <- add_par(doc2,
                paste0("Mean Root Mean Square Error of Item\u2013Total Correlation ",
                       "Estimators Across Test Lengths"),
                italic = TRUE, size = 12, pad_after = 6)
doc2 <- flextable::body_add_flextable(doc2, value = make_apa_ft_t56(rmse_tbl))
doc2 <- add_note(doc2,
                 paste0("k = number of items. Values are mean RMSE averaged over items and ",
                        R_sim, " simulation replications (N = ", N, " per replication). ",
                        "Lower values indicate closer agreement with the true item\u2013factor correlation. ",
                        "Abbreviations as in Table 5."))

doc2 <- officer::body_add_break(doc2)

# Figure 1: Bias plot
doc2 <- add_par(doc2, "Figure 1", bold = TRUE, size = 12, pad_after = 0)
doc2 <- add_par(doc2,
                paste0("Mean Bias of Six Item\u2013Total Correlation Estimators ",
                       "as a Function of Test Length"),
                italic = TRUE, size = 12, pad_after = 6)
doc2 <- officer::body_add_img(doc2, src = fig_bias_path, width = 6.5, height = 4.0)
doc2 <- add_note(doc2,
                 paste0("Mean bias across ", R_sim, " simulation replications (N = ", N,
                        " per replication). The horizontal reference line at zero denotes ",
                        "no bias. Lines below zero indicate underestimation; lines above zero ",
                        "indicate overestimation. Abbreviations as in Table 5."))

doc2 <- officer::body_add_break(doc2)

# Figure 2: RMSE plot
doc2 <- add_par(doc2, "Figure 2", bold = TRUE, size = 12, pad_after = 0)
doc2 <- add_par(doc2,
                paste0("Mean Root Mean Square Error of Six Item\u2013Total Correlation ",
                       "Estimators as a Function of Test Length"),
                italic = TRUE, size = 12, pad_after = 6)
doc2 <- officer::body_add_img(doc2, src = fig_rmse_path, width = 6.5, height = 4.0)
doc2 <- add_note(doc2,
                 paste0("Mean RMSE across ", R_sim, " simulation replications (N = ", N,
                        " per replication). Lower values indicate greater accuracy. ",
                        "Abbreviations as in Table 5."))

print(doc2, target = out_file_t56)
message(sprintf("Saved: %s", normalizePath(out_file_t56)))

message("\n=== ALL DONE ===")

################################################################################
# MANUSCRIPT NOTE:
# The column key in Section 3 (Results and Discussion) should be updated
# to reflect the new unified headers used across all six tables:
#
#   r-raw     -- uncorrected item-total correlation
#   Guilford  -- item-remainder correlation (Guilford, 1953, Eq. 1)
#   Cureton   -- corrected item-total, Cureton (1966, Eq. 3)
#   Henrysson -- corrected item-total, Henrysson (1963, Eq. 2)
#   CITC      -- proposed corrected item-total correlation (Eq. 11)
#   AIOSR     -- proposed estimator adjusted for overlap and scale
#                reliability (Eq. 12)
################################################################################