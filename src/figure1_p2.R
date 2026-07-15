# Figure 1 — Paper 2 (Gap 2)
# Secular trend in child stunting prevalence across five EDHS rounds (2000–2019)
# Two panels:
#   Left:  Observed weighted stunting prevalence (%) by survey year
#   Right: Adjusted survey-year odds ratios (DiD model, 2000 reference)
#          with crude ORs overlaid for comparison
# Vertical dashed reference line at 2013 (OWNP launch)
# Output: output/figures/figure1_p2.tiff (300 DPI, LZW) + figure1_p2.png (150 DPI)

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# ---- Load model outputs (no need to refit) ----------------------------------
full_tab   <- read.csv("output/tables/gap2_did_intensity_full.csv",
                       stringsAsFactors = FALSE)
crude_tab  <- read.csv("output/tables/crude_or_gap2_trend.csv",
                       stringsAsFactors = FALSE)
panel_data <- read.csv("data/data_clean/panel_5rounds_child.csv",
                       stringsAsFactors = FALSE, na.strings = c("NA",""))

# ---- Panel A: Weighted stunting prevalence by round -------------------------
intensity_tab <- panel_data %>%
  filter(survey_year == 2011, !is.na(edhs_improved_bin),
         !is.na(child_wt), child_wt > 0) %>%
  group_by(region) %>%
  summarise(intensity = sum((1 - edhs_improved_bin) * child_wt) / sum(child_wt),
            .groups = "drop")

model_df <- panel_data %>%
  left_join(intensity_tab %>% select(region, intensity), by = "region") %>%
  filter(!is.na(stunted_num), !is.na(child_wt), child_wt > 0,
         !is.na(edhs_improved_bin), !is.na(intensity))

# Weighted prevalence per round
prev_tab <- model_df %>%
  group_by(survey_year) %>%
  summarise(
    prev_pct = 100 * sum(stunted_num * child_wt, na.rm = TRUE) /
                     sum(child_wt, na.rm = TRUE),
    .groups = "drop"
  )

cat("=== Weighted stunting prevalence by round ===\n")
print(prev_tab)

# ---- Panel B: Survey-year ORs (crude vs adjusted, 2000 reference) -----------
# Adjusted: extract survey_year_f terms from DiD full model
adj_years <- full_tab %>%
  filter(grepl("^survey_year_f", term)) %>%
  mutate(
    survey_year = as.integer(gsub("survey_year_f", "", term)),
    type = "Adjusted"
  ) %>%
  select(survey_year, aOR, CI_low, CI_high, type)

# Add 2000 reference row (aOR = 1, CI = 1 to 1)
ref_adj <- data.frame(survey_year = 2000, aOR = 1, CI_low = 1, CI_high = 1, type = "Adjusted")
adj_years <- bind_rows(ref_adj, adj_years) %>% arrange(survey_year)

# Crude: from crude_or_gap2_trend.csv (2000 reference = Intercept)
crude_years <- crude_tab %>%
  filter(term != "(Intercept)") %>%
  mutate(
    survey_year = as.integer(gsub("survey_year_f", "", term)),
    type = "Crude (unadjusted)"
  ) %>%
  select(survey_year, aOR = cOR, CI_low, CI_high, type)

ref_crude <- data.frame(survey_year = 2000, aOR = 1, CI_low = 1, CI_high = 1,
                        type = "Crude (unadjusted)")
crude_years <- bind_rows(ref_crude, crude_years) %>% arrange(survey_year)

or_df <- bind_rows(adj_years, crude_years)
cat("\n=== Survey-year ORs ===\n")
print(or_df)

# ---- Build plots ------------------------------------------------------------
clrs <- c("Adjusted" = "#1B7EC2", "Crude (unadjusted)" = "#CC4125")
ltys  <- c("Adjusted" = "solid",  "Crude (unadjusted)" = "dashed")

# Panel A
pA <- ggplot(prev_tab, aes(x = survey_year, y = prev_pct)) +
  geom_vline(xintercept = 2013, linetype = "dotted", colour = "grey50", linewidth = 0.7) +
  geom_line(colour = "#333333", linewidth = 1.0) +
  geom_point(size = 3, colour = "#333333", fill = "white", shape = 21, stroke = 1.2) +
  geom_text(aes(label = sprintf("%.1f%%", prev_pct)),
            vjust = -1.0, size = 3.0, colour = "#333333") +
  annotate("text", x = 2013.3, y = 55, label = "OWNP\nlaunch\n2013",
           size = 2.8, hjust = 0, colour = "grey40", fontface = "italic") +
  scale_x_continuous(breaks = sort(unique(model_df$survey_year))) +
  scale_y_continuous(limits = c(30, 62), breaks = seq(30, 60, 10),
                     labels = function(y) paste0(y, "%")) +
  labs(x = "Survey year", y = "Child stunting prevalence (%)",
       title = "A  Observed weighted prevalence") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4))

# Panel B
pB <- ggplot(or_df, aes(x = survey_year, y = aOR, colour = type, linetype = type)) +
  geom_hline(yintercept = 1, colour = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 2013, linetype = "dotted", colour = "grey50", linewidth = 0.7) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = type),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  scale_linetype_manual(values = ltys) +
  scale_x_continuous(breaks = sort(unique(model_df$survey_year))) +
  scale_y_log10(limits = c(0.25, 2.2),
                breaks = c(0.3, 0.5, 0.75, 1.0, 1.5, 2.0),
                labels = c("0·30","0·50","0·75","1·00","1·50","2·00")) +
  labs(x = "Survey year", y = "Odds ratio (2000 reference; log scale)",
       colour = NULL, linetype = NULL, fill = NULL,
       title = "B  Survey-year odds ratios (reference: 2000)",
       caption = paste0(
         "Adjusted ORs from primary DiD model (region + survey-year fixed effects, wealth,\n",
         "maternal education, child age and sex). Crude ORs from bivariate svyglm (survey-year\n",
         "fixed effects only). Shaded bands: 95% CI. Dotted line: 2013 OWNP launch."
       )) +
  theme_classic(base_size = 10) +
  theme(legend.position    = "bottom",
        legend.text        = element_text(size = 9),
        plot.title         = element_text(size = 10, face = "bold"),
        plot.caption       = element_text(size = 7.5, colour = "grey40", hjust = 0),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4))

# ---- Combine with cowplot or patchwork (use grid if unavailable) -----------
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

if (has_patchwork) {
  library(patchwork)
  combined <- pA + pB + plot_layout(ncol = 2, widths = c(1, 1.2))
} else {
  library(gridExtra)
  combined <- arrangeGrob(pA, pB, ncol = 2, widths = c(1, 1.2))
}

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

if (has_patchwork) {
  ggsave("output/figures/figure1_p2.tiff", combined, width = 13, height = 6,
         dpi = 300, compression = "lzw")
  ggsave("output/figures/figure1_p2.png", combined, width = 13, height = 6, dpi = 150)
} else {
  tiff("output/figures/figure1_p2.tiff", width = 13, height = 6, units = "in",
       res = 300, compression = "lzw")
  grid::grid.draw(combined)
  dev.off()

  png("output/figures/figure1_p2.png", width = 13, height = 6, units = "in", res = 150)
  grid::grid.draw(combined)
  dev.off()
}

cat("Saved -> output/figures/figure1_p2.tiff\n")
cat("Saved -> output/figures/figure1_p2.png\n")
