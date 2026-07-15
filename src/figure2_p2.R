# Figure 2 — Paper 2 (Gap 2)
# Event-study forest plot: edhs_improved_bin × survey_year_f interaction coefficients
# (relative to 2011 reference), with 95% CIs
# This directly visualises the "no pre-trend or post-break" result
# alongside the primary DiD interaction terms as a summary panel
# Output: output/figures/figure2_p2.tiff (300 DPI, LZW) + figure2_p2.png (150 DPI)

suppressMessages({
  library(dplyr)
  library(ggplot2)
})

# ---- Load event-study coefficients (already saved from event_study_placebo_gap2.R)
es_tab <- read.csv("output/tables/gap2_event_study_coefs.csv",
                   stringsAsFactors = FALSE)

# ---- Load DiD interaction terms for summary panel
did_tab <- read.csv("output/tables/gap2_did_intensity_terms.csv",
                    stringsAsFactors = FALSE)

cat("=== Event-study coefficients ===\n")
print(es_tab)
cat("\n=== DiD interaction terms ===\n")
print(did_tab)

# ---- Panel A: Event-study plot -----------------------------------------------
# Exclude the reference row (2011, aOR=1)
es_plot <- es_tab %>%
  filter(!is.na(p_value) | survey_year == 2011) %>%
  mutate(
    label = ifelse(survey_year == 2011, "2011\n(ref)", as.character(survey_year)),
    is_ref = survey_year == 2011,
    pre_post = ifelse(survey_year < 2013, "Pre-OWNP", "Post-OWNP")
  )

# Add reference row explicitly (aOR=1, CI=1,1)
ref_row <- data.frame(
  term = "2011 (reference, omitted)", survey_year = 2011,
  beta = 0, ci_low_beta = 0, ci_high_beta = 0,
  aOR = 1, OR_CI_low = 1, OR_CI_high = 1, p_value = NA,
  label = "2011\n(ref)", is_ref = TRUE, pre_post = "Pre-OWNP"
)

es_plot <- es_tab %>%
  filter(survey_year != 2011) %>%
  mutate(
    label    = as.character(survey_year),
    is_ref   = FALSE,
    pre_post = ifelse(survey_year < 2013, "Pre-OWNP", "Post-OWNP")
  ) %>%
  bind_rows(ref_row) %>%
  arrange(survey_year)

clrs_es <- c("Pre-OWNP" = "#1B7EC2", "Post-OWNP" = "#CC4125")

es_plot <- es_plot %>% arrange(survey_year) %>% mutate(x_pos = row_number())
n_pre <- sum(es_plot$pre_post == "Pre-OWNP")
n_rounds_es <- nrow(es_plot)
post_xmin <- n_pre + 0.5
post_xmid <- (n_pre + 1 + n_rounds_es) / 2

pA <- ggplot(es_plot, aes(x = factor(survey_year), y = aOR,
                           colour = pre_post, shape = is_ref)) +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.5, linetype = "dashed") +
  geom_errorbar(aes(ymin = OR_CI_low, ymax = OR_CI_high),
                width = 0.15, linewidth = 0.7) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = clrs_es, name = "OWNP period") +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
  scale_y_log10(limits = c(0.25, 7),
                breaks = c(0.3, 0.5, 1, 2, 4),
                labels = c("0·30","0·50","1·00","2·00","4·00")) +
  annotate("rect", xmin = post_xmin, xmax = n_rounds_es + 0.5, ymin = 0.25, ymax = 7,
           alpha = 0.05, fill = "#CC4125") +
  annotate("text", x = post_xmid, y = 6.5,
           label = "Post-OWNP\nrounds", size = 2.8, colour = "#CC4125", fontface = "italic") +
  labs(
    x = "Survey year (interaction with improved water source use; 2011 reference)",
    y = "Odds ratio (log scale)",
    title = paste0(
      "A  Event-study: improved water source × survey-year interaction\n",
      "    (improved water–stunting association relative to 2011 baseline)"
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position    = "bottom",
    legend.text        = element_text(size = 9),
    plot.title         = element_text(size = 9.5, face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4)
  )

# ---- Recompute main-model N/PSU/strata for the panel B caption (avoids
# hardcoding a sample size that goes stale whenever a round is added) --------
panel_main <- read.csv("data/data_clean/panel_5rounds_child.csv",
                       stringsAsFactors = FALSE, na.strings = c("NA", ""))
intensity_main <- panel_main %>%
  filter(survey_year == 2011, !is.na(edhs_improved_bin), !is.na(child_wt), child_wt > 0) %>%
  group_by(region) %>%
  summarise(intensity = sum((1 - edhs_improved_bin) * child_wt) / sum(child_wt), .groups = "drop")
model_df_main <- panel_main %>%
  left_join(intensity_main, by = "region") %>%
  filter(!is.na(stunted_num), !is.na(child_wt), child_wt > 0,
         !is.na(edhs_improved_bin), !is.na(intensity))
n_main       <- nrow(model_df_main)
n_psu_main   <- length(unique(paste0(model_df_main$survey_year, "_", model_df_main$v001)))
n_strata_main <- length(unique(paste0(model_df_main$survey_year, "_", model_df_main$v023)))

# ---- Panel B: DiD interaction terms summary ---------------------------------
# Rename for readable plot labels
did_plot <- did_tab %>%
  mutate(term_label = case_when(
    term == "edhs_improved_bin"                      ~ "Improved water\n(main effect)",
    term == "edhs_improved_bin:post_ownp"             ~ "Water × Post-OWNP\n(two-way)",
    term == "edhs_improved_bin:intensity"             ~ "Water × Intensity\n(two-way)",
    term == "post_ownp:intensity"                    ~ "Post-OWNP × Intensity\n(two-way)",
    term == "edhs_improved_bin:post_ownp:intensity"  ~ "Water × Post-OWNP\n× Intensity\n(triple, primary)",
    TRUE ~ term
  )) %>%
  mutate(is_primary = grepl("triple", term_label))

did_plot <- did_plot %>%
  mutate(row_idx = seq_len(nrow(did_plot)))

pB <- ggplot(did_plot, aes(y = reorder(term_label, row_idx),
                            x = aOR, xmin = CI_low, xmax = CI_high,
                            colour = is_primary)) +
  geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.5, linetype = "dashed") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), height = 0.25,
                linewidth = 0.7, orientation = "y") +
  geom_point(size = 3.5) +
  geom_text(aes(label = sprintf("%.2f\n(%.2f–%.2f)", aOR, CI_low, CI_high)),
            nudge_x = 3, size = 2.6, hjust = 0) +
  scale_colour_manual(values = c("FALSE" = "grey40", "TRUE" = "#CC4125"), guide = "none") +
  scale_x_log10(limits = c(0.15, 20),
                breaks = c(0.2, 0.5, 1, 2, 5, 10),
                labels = c("0·20","0·50","1·00","2·00","5·00","10·0")) +
  labs(
    x = "Odds ratio (log scale)",
    y = NULL,
    title = "B  Generalised intensity DiD interaction terms\n    (primary triple interaction highlighted in red)",
    caption = paste0(
      "aOR = adjusted odds ratio; CI = 95% confidence interval.\n",
      "All interaction terms from primary DiD model (region + survey-year FEs, wealth, ",
      "education, child age and sex).\n",
      sprintf("N = %s children, %s round-qualified PSU-round observations, %s strata.",
              format(n_main, big.mark = ","), format(n_psu_main, big.mark = ","),
              format(n_strata_main, big.mark = ","))
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title   = element_text(size = 9.5, face = "bold"),
    plot.caption = element_text(size = 7.5, colour = "grey40", hjust = 0),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.4),
    axis.text.y  = element_text(size = 8.5)
  )

# ---- Combine ----------------------------------------------------------------
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

if (has_patchwork) {
  library(patchwork)
  combined <- pA + pB + plot_layout(ncol = 2, widths = c(1.1, 1))
} else {
  library(gridExtra)
  combined <- arrangeGrob(pA, pB, ncol = 2, widths = c(1.1, 1))
}

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

if (has_patchwork) {
  ggsave("output/figures/figure2_p2.tiff", combined, width = 14, height = 7,
         dpi = 300, compression = "lzw")
  ggsave("output/figures/figure2_p2.png", combined, width = 14, height = 7, dpi = 150)
} else {
  tiff("output/figures/figure2_p2.tiff", width = 14, height = 7, units = "in",
       res = 300, compression = "lzw")
  grid::grid.draw(combined)
  dev.off()
  png("output/figures/figure2_p2.png", width = 14, height = 7, units = "in", res = 150)
  grid::grid.draw(combined)
  dev.off()
}

cat("Saved -> output/figures/figure2_p2.tiff\n")
cat("Saved -> output/figures/figure2_p2.png\n")
