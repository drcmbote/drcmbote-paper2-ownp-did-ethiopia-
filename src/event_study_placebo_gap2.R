# Gap 2, pre-trends robustness: (1) event-study disaggregation of the
# intensity-DiD design around the 2011 reference round, and (2) a placebo
# DiD restricted to the two genuinely pre-OWNP rounds (2000, 2005), with
# 2005 relabelled as a fake "post" period.
#
# Both reuse the exact identification logic of src/did_gap2_intensity.R:
# baseline intensity is the 2011 (last pre-OWNP round) region-level,
# child_wt-weighted share of unimproved water source use; intensity is
# region-constant and therefore collinear with factor(region), so no
# standalone "intensity" main effect is included whenever region fixed
# effects are also in the model (identification comes only from terms where
# intensity is crossed with something that varies within region: round,
# or improved-water status). Round-qualified PSU/stratum IDs
# (survey_year_v001 / survey_year_v023) are used throughout because DHS
# cluster/stratum codes are not globally unique across rounds.

suppressMessages({
  library(survey)
  library(dplyr)
  library(ggplot2)
})

panel <- read.csv(
  "data/data_clean/panel_5rounds_child.csv",
  stringsAsFactors = FALSE, na.strings = c("NA", "")
)

intensity_2011 <- panel %>%
  filter(survey_year == 2011, !is.na(edhs_improved_bin), !is.na(child_wt), child_wt > 0) %>%
  group_by(region) %>%
  summarise(intensity = sum((1 - edhs_improved_bin) * child_wt) / sum(child_wt), .groups = "drop")

panel <- panel %>% left_join(intensity_2011, by = "region")

prep_covariates <- function(df) {
  df %>%
    mutate(
      region = factor(region),
      wealth_q = factor(wealth_q, levels = c("poorest", "poorer", "middle", "richer", "richest")),
      mother_education_f = factor(mother_education,
                                   levels = c("no_education", "primary", "secondary", "higher")),
      combined_psu = paste0(survey_year, "_", v001),
      combined_strata = paste0(survey_year, "_", v023)
    )
}

# =============================================================================
# [1] Event study: intensity x round-dummy, 2011 as reference period
# =============================================================================
cat("=== [1/2] Event study: intensity x survey-round, 2011 reference ===\n")

panel_es <- panel %>%
  prep_covariates() %>%
  mutate(
    survey_year_f = relevel(factor(survey_year), ref = "2011"),
    # Manual event-time dummies x intensity, 2011 deliberately omitted (=
    # reference, coefficient fixed at 0). R's default treatment contrasts
    # for "intensity:survey_year_f" would otherwise use FULL dummy coding
    # for survey_year_f in that interaction (because the standalone
    # "intensity" main effect is absent, dropped for its collinearity with
    # region) -- i.e. it would estimate a separate, non-relative intensity
    # slope for every round including 2011, rather than the intended
    # change-relative-to-2011 event-study coefficients. Building the
    # round-specific intensity columns by hand and leaving out int_2011
    # entirely is the standard way to force the 2011-omitted normalization.
    int_2000 = intensity * as.numeric(survey_year == 2000),
    int_2005 = intensity * as.numeric(survey_year == 2005),
    int_2016 = intensity * as.numeric(survey_year == 2016),
    int_2019 = intensity * as.numeric(survey_year == 2019),
    int_2024 = intensity * as.numeric(survey_year == 2024)
  )

model_es <- panel_es %>%
  filter(!is.na(stunted_num), !is.na(child_wt), child_wt > 0,
         !is.na(edhs_improved_bin), !is.na(intensity))

cat(sprintf("Event-study analytic sample: n=%d children, %d round-qualified PSUs\n",
            nrow(model_es), length(unique(model_es$combined_psu))))
print(table(model_es$survey_year))

des_es <- svydesign(ids = ~combined_psu, strata = ~combined_strata, weights = ~child_wt,
                     data = model_es, nest = TRUE)

fit_es <- svyglm(
  stunted_num ~ edhs_improved_bin + int_2000 + int_2005 + int_2016 + int_2019 + int_2024 +
    region + survey_year_f + wealth_q + child_age_months + child_sex_male_num + mother_education_f,
  design = des_es, family = quasibinomial(link = "logit")
)

s_es <- summary(fit_es)
beta <- coef(fit_es)
ci <- confint(fit_es)
p_val <- s_es$coefficients[, "Pr(>|t|)"]

es_terms <- c("int_2000", "int_2005", "int_2016", "int_2019", "int_2024")
es_tab <- data.frame(
  term = es_terms,
  survey_year = c(2000, 2005, 2016, 2019, 2024),
  beta = round(beta[es_terms], 4),
  ci_low_beta = round(ci[es_terms, 1], 4),
  ci_high_beta = round(ci[es_terms, 2], 4),
  aOR = round(exp(beta[es_terms]), 3),
  OR_CI_low = round(exp(ci[es_terms, 1]), 3),
  OR_CI_high = round(exp(ci[es_terms, 2]), 3),
  p_value = round(p_val[es_terms], 4)
)
ref_row <- data.frame(term = "2011 (reference, omitted)", survey_year = 2011,
                       beta = 0, ci_low_beta = 0, ci_high_beta = 0,
                       aOR = 1, OR_CI_low = 1, OR_CI_high = 1, p_value = NA)
es_tab <- rbind(es_tab, ref_row) %>% arrange(survey_year)

cat("\n--- Event-study coefficients: intensity x round (2011 = reference) ---\n")
print(es_tab, row.names = FALSE)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(es_tab, "output/tables/gap2_event_study_coefs.csv", row.names = FALSE)
cat("Saved -> output/tables/gap2_event_study_coefs.csv\n")

p_event <- ggplot(es_tab, aes(x = survey_year, y = beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 2011, linetype = "dotted", colour = "grey60") +
  geom_errorbar(aes(ymin = ci_low_beta, ymax = ci_high_beta), width = 0.6, colour = "steelblue") +
  geom_point(size = 2.8, colour = "steelblue") +
  scale_x_continuous(breaks = es_tab$survey_year) +
  labs(
    x = "EDHS survey round",
    y = "Intensity x round coefficient (log-odds, ref. = 2011)",
    title = "Gap 2: event-study test of pre-OWNP parallel trends",
    subtitle = "Regional baseline-intensity gradient in child stunting, by survey round"
  ) +
  theme_minimal(base_size = 13)

ggsave("output/figures/gap2_event_study_parallel.tiff", p_event,
       width = 8, height = 6, dpi = 300, compression = "lzw")
cat("Saved -> output/figures/gap2_event_study_parallel.tiff\n")

# =============================================================================
# [2] Placebo DiD: 2000 vs 2005 only, 2005 relabelled as a fake post period
# =============================================================================
cat("\n=== [2/2] Placebo DiD: 2000 (pre) vs 2005 (fake post) ===\n")

intensity_2000 <- panel %>%
  filter(survey_year == 2000, !is.na(edhs_improved_bin), !is.na(child_wt), child_wt > 0) %>%
  group_by(region) %>%
  summarise(intensity_2000 = sum((1 - edhs_improved_bin) * child_wt) / sum(child_wt), .groups = "drop")

panel_pl <- panel %>%
  filter(survey_year %in% c(2000, 2005)) %>%
  select(-intensity) %>%
  left_join(intensity_2000, by = "region") %>%
  rename(intensity = intensity_2000) %>%
  prep_covariates() %>%
  mutate(fake_post = ifelse(survey_year == 2005, 1, 0))

model_pl <- panel_pl %>%
  filter(!is.na(stunted_num), !is.na(child_wt), child_wt > 0,
         !is.na(edhs_improved_bin), !is.na(intensity))

cat(sprintf("Placebo analytic sample: n=%d children (2000+2005 only), %d round-qualified PSUs\n",
            nrow(model_pl), length(unique(model_pl$combined_psu))))
print(table(model_pl$survey_year))

des_pl <- svydesign(ids = ~combined_psu, strata = ~combined_strata, weights = ~child_wt,
                     data = model_pl, nest = TRUE)

fit_pl <- svyglm(
  stunted_num ~ edhs_improved_bin + edhs_improved_bin:fake_post + edhs_improved_bin:intensity +
    fake_post + fake_post:intensity + edhs_improved_bin:fake_post:intensity +
    region + wealth_q + child_age_months + child_sex_male_num + mother_education_f,
  design = des_pl, family = quasibinomial(link = "logit")
)

cat("\n========== Placebo DiD (2000 vs fake-post 2005) -- Child Stunting ==========\n")
print(summary(fit_pl))

or_tab <- exp(cbind(aOR = coef(fit_pl), confint(fit_pl)))
p_val_pl <- summary(fit_pl)$coefficients[, "Pr(>|t|)"]
tab_pl <- data.frame(
  term = rownames(or_tab), aOR = round(or_tab[, 1], 3),
  CI_low = round(or_tab[, 2], 3), CI_high = round(or_tab[, 3], 3),
  p_value = round(p_val_pl, 4)
)

placebo_terms <- tab_pl[tab_pl$term %in% c(
  "edhs_improved_bin", "edhs_improved_bin:fake_post", "edhs_improved_bin:intensity",
  "fake_post", "fake_post:intensity", "edhs_improved_bin:fake_post:intensity"
), ]
placebo_terms$n <- nrow(model_pl)

cat("\n--- Placebo DiD terms of interest ---\n")
print(placebo_terms, row.names = FALSE)

write.csv(placebo_terms, "output/tables/gap2_placebo_test.csv", row.names = FALSE)
cat("\nSaved -> output/tables/gap2_placebo_test.csv\n")

cat("\nDone.\n")
