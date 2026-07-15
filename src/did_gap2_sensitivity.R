# Gap 2 sensitivity analyses for the generalized intensity DiD
# (src/did_gap2_intensity.R is the main/reference model, already run and
# saved to output/tables/gap2_did_intensity_*.csv).
#
# Model A: drop 2019 (whose h11/diarrhea module is entirely missing, forcing
#   its exclusion from any diarrhea-adjusted model) and instead keep
#   recent_diarrhea_num as a covariate, using the rounds where it is fully
#   populated (2000/2005/2011/2016/2024-25 -- confirmed h11 is fielded again
#   in 2024-25, unlike the 2019 Interim round). Tests whether the main
#   model's null result was an artifact of dropping that covariate rather
#   than dropping 2019 itself.
#
# Model B: redefine treatment intensity using the EARLIEST round (2000)
#   region-level unimproved-water share instead of 2011 (the last pre-OWNP
#   round used in the main model), keeping all 5 rounds and the same
#   covariate set as the main model (no recent_diarrhea_num, for comparability
#   with the n=35,959 main-model sample). Tests sensitivity of the DiD
#   estimate to the choice of baseline year for the intensity measure.
#
# Both reuse the same round-qualified PSU/stratum IDs (survey_year_v001 /
# survey_year_v023) as the main model, since DHS cluster/stratum codes are
# not globally unique across rounds.

suppressMessages({
  library(survey)
  library(dplyr)
})

panel_raw <- read.csv(
  "data/data_clean/panel_5rounds_child.csv",
  stringsAsFactors = FALSE, na.strings = c("NA", "")
)

build_intensity <- function(panel, baseline_year) {
  panel %>%
    filter(survey_year == baseline_year, !is.na(edhs_improved_bin), !is.na(child_wt), child_wt > 0) %>%
    group_by(region) %>%
    summarise(intensity = sum((1 - edhs_improved_bin) * child_wt) / sum(child_wt), .groups = "drop")
}

prep_covariates <- function(panel) {
  panel %>%
    mutate(
      region = factor(region),
      survey_year_f = factor(survey_year),
      wealth_q = factor(wealth_q, levels = c("poorest", "poorer", "middle", "richer", "richest")),
      mother_education_f = factor(mother_education,
                                   levels = c("no_education", "primary", "secondary", "higher")),
      combined_psu = paste0(survey_year, "_", v001),
      combined_strata = paste0(survey_year, "_", v023)
    )
}

extract_did_terms <- function(fit, model_label, n, n_rounds) {
  or_tab <- exp(cbind(aOR = coef(fit), confint(fit)))
  p_val <- summary(fit)$coefficients[, "Pr(>|t|)"]
  tab <- data.frame(
    term = rownames(or_tab), aOR = round(or_tab[, 1], 3),
    CI_low = round(or_tab[, 2], 3), CI_high = round(or_tab[, 3], 3),
    p_value = round(p_val, 4)
  )
  did <- tab[tab$term %in% c("edhs_improved_bin", "edhs_improved_bin:post_ownp",
                              "edhs_improved_bin:intensity", "post_ownp:intensity",
                              "edhs_improved_bin:post_ownp:intensity"), ]
  did$model <- model_label
  did$n <- n
  did$n_rounds <- n_rounds
  did[, c("model", "n", "n_rounds", "term", "aOR", "CI_low", "CI_high", "p_value")]
}

# ============================================================
# Model A: 4 rounds (drop 2019), diarrhea-adjusted, 2011 baseline intensity
# ============================================================
cat("=== Model A: rounds excluding 2019, recent_diarrhea_num included, 2011 baseline intensity ===\n")

intensity_2011 <- build_intensity(panel_raw, 2011)
panelA <- panel_raw %>%
  filter(survey_year != 2019) %>%
  left_join(intensity_2011, by = "region") %>%
  prep_covariates()

model_dfA <- panelA %>%
  filter(!is.na(stunted_num), !is.na(child_wt), child_wt > 0, !is.na(recent_diarrhea_num),
         !is.na(edhs_improved_bin), !is.na(intensity))

cat(sprintf("Model A analytic sample: n=%d, rounds=%s\n", nrow(model_dfA),
            paste(sort(unique(model_dfA$survey_year)), collapse = ",")))

desA <- svydesign(ids = ~combined_psu, strata = ~combined_strata, weights = ~child_wt,
                   data = model_dfA, nest = TRUE)

fitA <- svyglm(
  stunted_num ~ edhs_improved_bin + edhs_improved_bin:post_ownp + edhs_improved_bin:intensity +
    post_ownp:intensity + edhs_improved_bin:post_ownp:intensity +
    region + survey_year_f + wealth_q + child_age_months + child_sex_male_num +
    recent_diarrhea_num + mother_education_f,
  design = desA, family = quasibinomial(link = "logit")
)

roundsA <- paste(sort(unique(model_dfA$survey_year)), collapse = ",")
didA <- extract_did_terms(fitA, sprintf("Model A: %d rounds (no 2019) + diarrhea-adjusted, 2011 intensity",
                                         length(unique(model_dfA$survey_year))),
                           nrow(model_dfA), roundsA)
cat("\n--- Model A DiD terms ---\n")
print(didA, row.names = FALSE)

# ============================================================
# Model B: 5 rounds, 2000-baseline intensity, same covariates as main model
# ============================================================
cat("\n=== Model B: all 5 rounds, intensity redefined from 2000 baseline ===\n")

intensity_2000 <- build_intensity(panel_raw, 2000)
panelB <- panel_raw %>%
  left_join(intensity_2000, by = "region") %>%
  prep_covariates()

model_dfB <- panelB %>%
  filter(!is.na(stunted_num), !is.na(child_wt), child_wt > 0,
         !is.na(edhs_improved_bin), !is.na(intensity))

cat(sprintf("Model B analytic sample: n=%d, rounds=%s\n", nrow(model_dfB),
            paste(sort(unique(model_dfB$survey_year)), collapse = ",")))

desB <- svydesign(ids = ~combined_psu, strata = ~combined_strata, weights = ~child_wt,
                   data = model_dfB, nest = TRUE)

fitB <- svyglm(
  stunted_num ~ edhs_improved_bin + edhs_improved_bin:post_ownp + edhs_improved_bin:intensity +
    post_ownp:intensity + edhs_improved_bin:post_ownp:intensity +
    region + survey_year_f + wealth_q + child_age_months + child_sex_male_num +
    mother_education_f,
  design = desB, family = quasibinomial(link = "logit")
)

roundsB <- paste(sort(unique(model_dfB$survey_year)), collapse = ",")
didB <- extract_did_terms(fitB, sprintf("Model B: %d rounds, 2000-baseline intensity",
                                         length(unique(model_dfB$survey_year))),
                           nrow(model_dfB), roundsB)
cat("\n--- Model B DiD terms ---\n")
print(didB, row.names = FALSE)

# ============================================================
# Comparison table: Main + Model A + Model B
# ============================================================
main_df <- panel_raw %>%
  left_join(intensity_2011, by = "region") %>%
  filter(!is.na(stunted_num), !is.na(child_wt), child_wt > 0,
         !is.na(edhs_improved_bin), !is.na(intensity))
main_rounds <- paste(sort(unique(main_df$survey_year)), collapse = ",")

main_tab <- read.csv("output/tables/gap2_did_intensity_terms.csv", stringsAsFactors = FALSE)
main_tab$model <- sprintf("Main: %d rounds, 2011-baseline intensity", length(unique(main_df$survey_year)))
main_tab$n <- nrow(main_df)
main_tab$n_rounds <- main_rounds
main_tab <- main_tab[, c("model", "n", "n_rounds", "term", "aOR", "CI_low", "CI_high", "p_value")]

comparison <- rbind(main_tab, didA, didB)

cat("\n========== SIDE-BY-SIDE COMPARISON: Main vs Model A vs Model B ==========\n")
print(comparison, row.names = FALSE)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, "output/tables/gap2_sensitivity_comparison.csv", row.names = FALSE)
cat("\nSaved -> output/tables/gap2_sensitivity_comparison.csv\n")
