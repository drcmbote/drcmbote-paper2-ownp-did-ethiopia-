# STROBE Item 16(a): Unadjusted (crude) OR estimates for Paper 2 (Gap 2)
# Model: svyglm with complex-sampling weights
# Region and survey-year FEs retained (identification requirement, not confounders)
# No individual-level covariates (wealth, age, sex, maternal education omitted)
# DiD interaction terms preserved for structural comparability with adjusted model

suppressMessages({
  library(survey)
  library(dplyr)
})

panel <- read.csv(
  "data/data_clean/panel_5rounds_child.csv",
  stringsAsFactors = FALSE, na.strings = c("NA", "")
)

# Bartik intensity (2011 baseline) — identical to main did_gap2_intensity.R
intensity_tab <- panel %>%
  filter(survey_year == 2011, !is.na(edhs_improved_bin),
         !is.na(child_wt), child_wt > 0) %>%
  group_by(region) %>%
  summarise(intensity = sum((1 - edhs_improved_bin) * child_wt) / sum(child_wt))

panel <- panel %>%
  left_join(intensity_tab %>% select(region, intensity), by = "region") %>%
  mutate(
    region        = factor(region),
    survey_year_f = factor(survey_year),
    post_ownp     = ifelse(survey_year >= 2016, 1, 0),
    combined_psu    = paste0(survey_year, "_", v001),
    combined_strata = paste0(survey_year, "_", v023)
  )

model_df <- panel %>%
  filter(!is.na(stunted_num), !is.na(child_wt), child_wt > 0,
         !is.na(edhs_improved_bin), !is.na(intensity))

cat(sprintf("Crude DiD analytic sample: n=%d children\n", nrow(model_df)))

des <- svydesign(
  ids     = ~combined_psu,
  strata  = ~combined_strata,
  weights = ~child_wt,
  data    = model_df,
  nest    = TRUE
)

# Crude model: region + year FEs + DiD interaction terms; no individual covariates
fit_crude <- svyglm(
  stunted_num ~ edhs_improved_bin +
    edhs_improved_bin:post_ownp +
    edhs_improved_bin:intensity +
    post_ownp:intensity +
    edhs_improved_bin:post_ownp:intensity +
    region + survey_year_f,
  design = des, family = quasibinomial(link = "logit")
)

b   <- coef(fit_crude)
se  <- sqrt(diag(vcov(fit_crude)))
z   <- b / se
p   <- 2 * pt(-abs(z), df = fit_crude$df.residual)

crude_results <- data.frame(
  term    = names(b),
  cOR     = round(exp(b),              3),
  CI_low  = round(exp(b - 1.96 * se), 3),
  CI_high = round(exp(b + 1.96 * se), 3),
  p_value = round(p,                  4)
)

cat("\n========== CRUDE svyglm — Child Stunting (region + year FEs, no individual covariates) ==========\n")
print(crude_results, row.names = FALSE)

# ---- DiD terms of interest (primary table rows) ----------------------------
did_terms <- c("edhs_improved_bin",
               "edhs_improved_bin:post_ownp",
               "edhs_improved_bin:intensity",
               "post_ownp:intensity",
               "edhs_improved_bin:post_ownp:intensity")
cat("\n--- Crude DiD interaction terms ---\n")
print(crude_results[crude_results$term %in% did_terms, ], row.names = FALSE)

# ---- Survey-year trend (crude) ----------------------------------------------
# Simple bivariate stunting ~ survey_year_f with complex sampling (no FEs or covariates)
fit_trend_crude <- svyglm(
  stunted_num ~ survey_year_f,
  design = des, family = quasibinomial(link = "logit")
)

b2  <- coef(fit_trend_crude)
se2 <- sqrt(diag(vcov(fit_trend_crude)))
z2  <- b2 / se2
p2  <- 2 * pt(-abs(z2), df = fit_trend_crude$df.residual)

trend_crude <- data.frame(
  term    = names(b2),
  cOR     = round(exp(b2),               3),
  CI_low  = round(exp(b2 - 1.96 * se2), 3),
  CI_high = round(exp(b2 + 1.96 * se2), 3),
  p_value = round(p2,                   4)
)

cat("\n--- Crude survey-year ORs (bivariate, 2000 reference) ---\n")
print(trend_crude, row.names = FALSE)

# ---- Save -------------------------------------------------------------------
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(crude_results,  "output/tables/crude_or_gap2_did.csv",   row.names = FALSE)
write.csv(trend_crude,    "output/tables/crude_or_gap2_trend.csv", row.names = FALSE)
cat("\nSaved -> output/tables/crude_or_gap2_did.csv\n")
cat("Saved -> output/tables/crude_or_gap2_trend.csv\n")
