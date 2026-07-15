# Gap 2: generalized intensity difference-in-differences for the pre/post
# OWNP (national, simultaneous 2013 launch) design on the 5-round pooled
# EDHS child panel (2000/2005/2011/2016/2019).
#
# Treatment intensity: each region's 2011 (last pre-OWNP round) weighted
# share of households using an UNIMPROVED water source (child_wt-weighted).
# This is a standard Bartik/shift-share intensity measure -- regions with
# more unimproved-water exposure at baseline had more "room" for a national
# WASH program to matter. Intensity is constant within region across all 5
# survey rounds (it is mapped from the single 2011 cross-section onto every
# round by region).
#
# Identification note on collinearity: intensity is a deterministic,
# region-constant value, so it is perfectly collinear with region fixed
# effects (factor(region)); likewise post_ownp is a deterministic function
# of survey_year, perfectly collinear with year fixed effects
# (factor(survey_year)). Both standalone main effects are therefore dropped
# from the formula -- only the interaction terms involving edhs_improved_bin
# (a child/household-level variable that varies WITHIN region-year cells)
# survive as identified coefficients. This is the standard generalized/
# intensity-DiD specification (identification comes from how the post-vs-pre
# *change* in the improved-water/stunting relationship differs by baseline
# regional intensity, conditional on region and year fixed effects).
#
# Complex-sampling note: DHS cluster (v001) and stratum (v023) codes are
# NOT globally unique across survey rounds (cluster "1" in 2000 is a
# different physical PSU than cluster "1" in 2016), so the PSU/stratum IDs
# used in svydesign are round-qualified: combined_psu = survey_year_v001,
# combined_strata = survey_year_v023.
#
# Covariate note: recent_diarrhea_num (h11) is verified to be ENTIRELY
# missing in the 2019 Interim DHS KR file (the health-recall module was not
# fielded/released in that interim round) -- requiring it as a complete-case
# covariate would silently drop all 5,753 children from 2019, removing half
# of the post-OWNP arm. It is therefore excluded from this pooled panel
# model (it remains available and was already used in the single-round Gap 1
# models). child_age_months, child_sex_male_num and mother_education are
# confirmed 100% complete in all 5 rounds and are retained.

suppressMessages({
  library(survey)
  library(dplyr)
})

panel <- read.csv(
  "data/data_clean/panel_5rounds_child.csv",
  stringsAsFactors = FALSE, na.strings = c("NA", "")
)

# ---- 1. Treatment intensity: 2011 region-level unimproved-water share -----
intensity_tab <- panel %>%
  filter(survey_year == 2011, !is.na(edhs_improved_bin), !is.na(child_wt), child_wt > 0) %>%
  group_by(region) %>%
  summarise(
    intensity = sum((1 - edhs_improved_bin) * child_wt) / sum(child_wt),
    n_2011 = n()
  )

cat("=== Treatment intensity (2011 weighted unimproved-water share by region) ===\n")
print(as.data.frame(intensity_tab))

panel <- panel %>% left_join(intensity_tab %>% select(region, intensity), by = "region")

# ---- 2. Covariate construction --------------------------------------------
panel <- panel %>%
  mutate(
    region = factor(region),
    survey_year_f = factor(survey_year),
    wealth_q = factor(wealth_q, levels = c("poorest", "poorer", "middle", "richer", "richest")),
    mother_education_f = factor(mother_education,
                                 levels = c("no_education", "primary", "secondary", "higher")),
    combined_psu = paste0(survey_year, "_", v001),
    combined_strata = paste0(survey_year, "_", v023)
  )

model_df <- panel %>%
  filter(!is.na(stunted_num), !is.na(child_wt), child_wt > 0,
         !is.na(edhs_improved_bin), !is.na(intensity))

cat(sprintf("\nDiD analytic sample: n=%d children, %d round-qualified PSUs, %d round-qualified strata\n",
            nrow(model_df), length(unique(model_df$combined_psu)), length(unique(model_df$combined_strata))))
cat("\nRows per round in analytic sample:\n")
print(table(model_df$survey_year))

des <- svydesign(
  ids = ~combined_psu, strata = ~combined_strata, weights = ~child_wt,
  data = model_df, nest = TRUE
)

fit <- svyglm(
  stunted_num ~ edhs_improved_bin + edhs_improved_bin:post_ownp + edhs_improved_bin:intensity +
    post_ownp:intensity + edhs_improved_bin:post_ownp:intensity +
    region + survey_year_f + wealth_q + child_age_months + child_sex_male_num +
    mother_education_f,
  design = des, family = quasibinomial(link = "logit")
)

cat("\n========== Generalized Intensity DiD -- Child Stunting (HAZ < -2) ==========\n")
print(summary(fit))

or_tab <- exp(cbind(aOR = coef(fit), confint(fit)))
p_val <- summary(fit)$coefficients[, "Pr(>|t|)"]
tab <- data.frame(
  term = rownames(or_tab), aOR = round(or_tab[, 1], 3),
  CI_low = round(or_tab[, 2], 3), CI_high = round(or_tab[, 3], 3),
  p_value = round(p_val, 4)
)

cat("\n--- aOR (95% CI), p-value: full model ---\n")
print(tab, row.names = FALSE)

cat("\n--- aOR (95% CI), p-value: DiD terms of interest only ---\n")
did_terms <- tab[grepl("edhs_improved_bin|post_ownp|intensity", tab$term, fixed = FALSE) &
                   !grepl("^region|^survey_year_f", tab$term), ]
print(did_terms, row.names = FALSE)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(tab, "output/tables/gap2_did_intensity_full.csv", row.names = FALSE)
write.csv(did_terms, "output/tables/gap2_did_intensity_terms.csv", row.names = FALSE)
cat("\nSaved -> output/tables/gap2_did_intensity_full.csv\n")
cat("Saved -> output/tables/gap2_did_intensity_terms.csv\n")
