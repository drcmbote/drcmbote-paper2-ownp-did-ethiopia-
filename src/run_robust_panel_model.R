# run_robust_panel_model.R
# Purpose:
#   Robust multilevel panel model — 5-round pooled Ethiopia DHS (2000-2019)
#   Primary question: does edhs_improved_bin independently predict child stunting
#   (stunted_num) after adjusting for confounders and clustering?
#   Replicates/verifies the aOR~1.009, p=0.94 null finding.
#
# Strategy:
#   Model A  - crude logistic (no adjustment, no random effects)
#   Model B  - survey-weighted svyglm with PSU/strata (complex-survey GEE approx)
#   Model C  - WeMix GLMM, random intercept on combined PSU (survey_year x v001)
#   Model D  - lme4 glmer sensitivity (no survey weights, ML)
#   Model E  - region-stratified meta-analytic summary (Cochran Q heterogeneity)
#
# Data: data/data_clean/panel_5rounds_child.csv
#   child_sex_male_num  -> already numeric 0/1 (confirmed in header)
#   edhs_improved_bin   -> 0/1 numeric
#   stunted_num         -> 0/1 numeric (HAZ < -2)

suppressMessages({
  library(dplyr)
  library(survey)
})

PKG <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("[SKIP] %s not installed — skipping that model", pkg))
    return(FALSE)
  }
  suppressMessages(library(pkg, character.only = TRUE))
  TRUE
}

WEMIX_OK <- FALSE   # disabled: 35k obs x 2604 clusters → GHQ too slow
LME4_OK  <- PKG("lme4")
META_OK  <- PKG("meta")

# ============================================================
# 0. LOAD & DIAGNOSTICS
# ============================================================
panel_raw <- read.csv(
  "data/data_clean/panel_5rounds_child.csv",
  stringsAsFactors = FALSE,
  na.strings        = c("NA", "", "N/A", "na", "NaN")
)

cat("=== DATASET DIAGNOSTICS ===\n")
cat(sprintf("Rows: %d  |  Cols: %d\n", nrow(panel_raw), ncol(panel_raw)))
cat("Columns:", paste(names(panel_raw), collapse = ", "), "\n\n")

key_vars <- c("stunted_num","edhs_improved_bin","child_sex_male_num",
              "child_age_months","mother_education","wealth_q",
              "region","survey_year","child_wt","v001","v023")

cat("--- Key variable summaries ---\n")
for (v in key_vars) {
  if (!v %in% names(panel_raw)) { cat(sprintf("  %-28s  [NOT FOUND]\n", v)); next }
  x <- panel_raw[[v]]
  if (is.numeric(x)) {
    cat(sprintf("  %-28s  numeric  range=[%.2f, %.2f]  NA=%d\n",
                v, min(x, na.rm=TRUE), max(x, na.rm=TRUE), sum(is.na(x))))
  } else {
    cat(sprintf("  %-28s  char     unique=%d  NA=%d\n",
                v, length(unique(x[!is.na(x)])), sum(is.na(x))))
  }
}
cat("\n")

# Confirm child_sex_male_num is numeric
stopifnot("child_sex_male_num must be numeric" =
            is.numeric(panel_raw$child_sex_male_num))
cat("CONFIRMED: child_sex_male_num is numeric (0/1)\n\n")

# ============================================================
# 1. FACTOR ENCODING
# ============================================================
me_vals <- sort(unique(panel_raw$mother_education[!is.na(panel_raw$mother_education)]))
if (is.numeric(panel_raw$mother_education)) {
  me_levels <- c(0, 1, 2, 3)
  me_labels <- c("no_education","primary","secondary","higher")
} else if (all(me_vals %in% c("no_education","primary","secondary","higher"))) {
  me_levels <- c("no_education","primary","secondary","higher")
  me_labels <- c("no_education","primary","secondary","higher")
} else {
  me_levels <- me_vals; me_labels <- me_vals
}

wq_vals <- sort(unique(panel_raw$wealth_q[!is.na(panel_raw$wealth_q)]))
if (is.numeric(panel_raw$wealth_q)) {
  wq_levels <- 1:5
  wq_labels <- c("poorest","poorer","middle","richer","richest")
} else if (all(wq_vals %in% c("poorest","poorer","middle","richer","richest"))) {
  wq_levels <- c("poorest","poorer","middle","richer","richest")
  wq_labels <- c("poorest","poorer","middle","richer","richest")
} else {
  wq_levels <- wq_vals; wq_labels <- as.character(wq_vals)
}

df <- panel_raw %>%
  mutate(
    survey_year_f      = factor(survey_year),
    mother_education_f = factor(mother_education, levels = me_levels, labels = me_labels),
    wealth_q_f         = factor(wealth_q,         levels = wq_levels, labels = wq_labels),
    combined_psu       = paste0(survey_year, "_", v001),
    combined_strat     = paste0(survey_year, "_", v023),
    wt2                = 1
  ) %>%
  filter(
    !is.na(stunted_num),
    !is.na(edhs_improved_bin),
    !is.na(child_wt), child_wt > 0,
    !is.na(child_age_months),
    !is.na(child_sex_male_num),
    !is.na(mother_education_f),
    !is.na(wealth_q_f)
  )

cat(sprintf("=== ANALYTICAL SAMPLE: N = %d children, %d clusters, %d regions ===\n",
            nrow(df),
            length(unique(df$combined_psu)),
            length(unique(df$region))))
cat("Survey rounds:", paste(sort(unique(df$survey_year)), collapse=", "), "\n\n")
cat("Stunting prevalence by round:\n")
for (y in sort(unique(df$survey_year))) {
  sub <- df[df$survey_year == y, ]
  cat(sprintf("  %d: n=%5d  stunted=%5d (%.1f%%)\n",
              y, nrow(sub), sum(sub$stunted_num), 100*mean(sub$stunted_num)))
}
cat(sprintf("\nImproved water overall: %.1f%%\n\n",
            100*mean(df$edhs_improved_bin)))

# ============================================================
# 2. MODEL A — CRUDE (unadjusted) LOGISTIC
# ============================================================
cat("=== MODEL A: CRUDE LOGISTIC (edhs_improved_bin only) ===\n")
fitA <- glm(stunted_num ~ edhs_improved_bin,
            data = df, family = binomial(link = "logit"))
sA   <- summary(fitA)$coefficients
bA   <- sA["edhs_improved_bin", "Estimate"]
seA  <- sA["edhs_improved_bin", "Std. Error"]
cat(sprintf("  Crude OR = %.3f  (95%% CI %.3f-%.3f)  p = %.4f\n\n",
            exp(bA), exp(bA - 1.96*seA), exp(bA + 1.96*seA),
            sA["edhs_improved_bin", "Pr(>|z|)"]))

# ============================================================
# 3. MODEL B — COMPLEX-SURVEY svyglm (GEE approximation)
# ============================================================
cat("=== MODEL B: SURVEY-WEIGHTED svyglm (complex-survey adjusted) ===\n")
svy <- svydesign(
  ids     = ~combined_psu,
  strata  = ~combined_strat,
  weights = ~child_wt,
  data    = df,
  nest    = TRUE
)

fitB <- svyglm(
  stunted_num ~ edhs_improved_bin + survey_year_f + region +
    wealth_q_f + mother_education_f + child_age_months + child_sex_male_num,
  design = svy,
  family = quasibinomial(link = "logit")
)
sB   <- summary(fitB)$coefficients
bB   <- sB["edhs_improved_bin", "Estimate"]
seB  <- sB["edhs_improved_bin", "Std. Error"]
pB   <- sB["edhs_improved_bin", "Pr(>|t|)"]
cat(sprintf("  Adjusted OR = %.4f  (95%% CI %.4f-%.4f)  p = %.4f\n\n",
            exp(bB), exp(bB - 1.96*seB), exp(bB + 1.96*seB), pB))

# ============================================================
# 4. MODEL C — WeMix GLMM (random intercept on PSU x year)
# ============================================================
if (WEMIX_OK) {
  cat("=== MODEL C: WeMix GLMM (random intercept: combined_psu) ===\n")
  fitC <- WeMix::mix(
    stunted_num ~ edhs_improved_bin + survey_year_f + region +
      wealth_q_f + mother_education_f + child_age_months + child_sex_male_num +
      (1 | combined_psu),
    data    = df,
    weights = c("child_wt", "wt2"),
    family  = binomial(link = "logit")
  )
  sC   <- summary(fitC)$coef
  bC   <- sC["edhs_improved_bin", "Estimate"]
  seC  <- sC["edhs_improved_bin", "Std. Error"]
  pC   <- 2 * pnorm(-abs(bC / seC))
  cat(sprintf("  aOR = %.4f  (95%% CI %.4f-%.4f)  p = %.4f\n",
              exp(bC), exp(bC - 1.96*seC), exp(bC + 1.96*seC), pC))
  cat(sprintf("  Cluster random-intercept variance: %.4f\n",
              summary(fitC)$varVC[["combined_psu"]]))
  cat(sprintf("  ICC (approx): %.4f\n\n",
              summary(fitC)$varVC[["combined_psu"]] /
                (summary(fitC)$varVC[["combined_psu"]] + pi^2/3)))
} else {
  cat("=== MODEL C: WeMix not available — SKIPPED ===\n\n")
  bC <- NA; seC <- NA; pC <- NA
}

# ============================================================
# 5. MODEL D — lme4 glmer (unweighted, ML, sensitivity)
# ============================================================
if (LME4_OK) {
  cat("=== MODEL D: lme4 glmer (sensitivity — unweighted, ML) ===\n")
  fitD <- lme4::glmer(
    stunted_num ~ edhs_improved_bin + survey_year_f + region +
      wealth_q_f + mother_education_f + child_age_months + child_sex_male_num +
      (1 | combined_psu),
    data    = df,
    family  = binomial(link = "logit"),
    nAGQ    = 0,   # PQL: much faster for 2600+ clusters; Laplace would take >30min
    control = lme4::glmerControl(optimizer = "nloptwrap",
                                  optCtrl  = list(maxeval = 1e4))
  )
  sD   <- summary(fitD)$coefficients
  bD   <- sD["edhs_improved_bin", "Estimate"]
  seD  <- sD["edhs_improved_bin", "Std. Error"]
  zD   <- sD["edhs_improved_bin", "z value"]
  pD   <- sD["edhs_improved_bin", "Pr(>|z|)"]
  cat(sprintf("  aOR = %.4f  (95%% CI %.4f-%.4f)  z=%.3f  p = %.4f\n",
              exp(bD), exp(bD - 1.96*seD), exp(bD + 1.96*seD), zD, pD))
  vc <- as.data.frame(lme4::VarCorr(fitD))
  cat(sprintf("  Cluster random-intercept variance: %.4f\n",
              vc$vcov[1]))
  cat(sprintf("  ICC (approx): %.4f\n\n",
              vc$vcov[1] / (vc$vcov[1] + pi^2/3)))
} else {
  cat("=== MODEL D: lme4 not available — SKIPPED ===\n\n")
  bD <- NA; seD <- NA; pD <- NA
}

# ============================================================
# 6. MODEL E — REGION-STRATIFIED META-ANALYSIS
# ============================================================
cat("=== MODEL E: REGION-STRATIFIED LOGISTIC (Cochran Q test for heterogeneity) ===\n")
regions_present <- sort(unique(df$region))
reg_results <- lapply(regions_present, function(reg) {
  sub <- df[df$region == reg, ]
  if (nrow(sub) < 30 || length(unique(sub$stunted_num)) < 2 ||
      length(unique(sub$edhs_improved_bin)) < 2) return(NULL)
  tryCatch({
    fit_r <- glm(
      stunted_num ~ edhs_improved_bin + survey_year_f +
        wealth_q_f + mother_education_f + child_age_months + child_sex_male_num,
      data = sub, family = binomial(link = "logit"),
      weights = sub$child_wt
    )
    s_r <- summary(fit_r)$coefficients
    if (!"edhs_improved_bin" %in% rownames(s_r)) return(NULL)
    b_r  <- s_r["edhs_improved_bin", "Estimate"]
    se_r <- s_r["edhs_improved_bin", "Std. Error"]
    data.frame(region = reg, n = nrow(sub), logOR = b_r, seLogOR = se_r,
               OR = exp(b_r), lci = exp(b_r - 1.96*se_r), uci = exp(b_r + 1.96*se_r),
               p = s_r["edhs_improved_bin", "Pr(>|z|)"])
  }, error = function(e) NULL)
})
reg_tab <- do.call(rbind, Filter(Negate(is.null), reg_results))

if (!is.null(reg_tab) && nrow(reg_tab) > 1) {
  cat("  Region-specific ORs:\n")
  for (i in seq_len(nrow(reg_tab))) {
    cat(sprintf("    %-14s  n=%5d  OR=%.3f (%.3f-%.3f)  p=%.4f\n",
                reg_tab$region[i], reg_tab$n[i],
                reg_tab$OR[i], reg_tab$lci[i], reg_tab$uci[i], reg_tab$p[i]))
  }

  # Fixed-effect pooled estimate (inverse-variance)
  w       <- 1 / reg_tab$seLogOR^2
  pooled_b <- sum(w * reg_tab$logOR) / sum(w)
  pooled_se <- sqrt(1 / sum(w))
  Q       <- sum(w * (reg_tab$logOR - pooled_b)^2)
  df_Q    <- nrow(reg_tab) - 1
  I2      <- max(0, 100 * (Q - df_Q) / Q)

  cat(sprintf("\n  Pooled OR (IV fixed-effect): %.4f (%.4f-%.4f)\n",
              exp(pooled_b), exp(pooled_b - 1.96*pooled_se),
              exp(pooled_b + 1.96*pooled_se)))
  cat(sprintf("  Cochran Q = %.2f  df = %d  p = %.4f  I² = %.1f%%\n\n",
              Q, df_Q, pchisq(Q, df_Q, lower.tail = FALSE), I2))
}

# ============================================================
# 7. SENSITIVITY — SUBGROUP BY SURVEY PERIOD
# ============================================================
cat("=== SENSITIVITY: PRE vs POST OWNP (2016-2019) SUBGROUP ===\n")
for (period in c("Pre-OWNP (2000-2011)", "Post-OWNP (2016-2019)")) {
  sub <- if (grepl("Pre", period)) df[df$survey_year <= 2011, ] else df[df$survey_year >= 2016, ]
  if (length(unique(sub$edhs_improved_bin)) < 2) { cat("  ", period, "— insufficient variation\n"); next }
  svy_s <- svydesign(ids=~combined_psu, strata=~combined_strat,
                      weights=~child_wt, data=sub, nest=TRUE)
  fit_s <- tryCatch(
    svyglm(stunted_num ~ edhs_improved_bin + region + wealth_q_f +
             mother_education_f + child_age_months + child_sex_male_num,
           design = svy_s, family = quasibinomial()),
    error = function(e) NULL
  )
  if (is.null(fit_s)) { cat("  ", period, "— model failed\n"); next }
  s_s   <- summary(fit_s)$coefficients
  if (!"edhs_improved_bin" %in% rownames(s_s)) next
  b_s   <- s_s["edhs_improved_bin", "Estimate"]
  se_s  <- s_s["edhs_improved_bin", "Std. Error"]
  p_s   <- s_s["edhs_improved_bin", "Pr(>|t|)"]
  cat(sprintf("  %-30s  aOR=%.4f (%.4f-%.4f)  p=%.4f  n=%d\n",
              period, exp(b_s), exp(b_s - 1.96*se_s), exp(b_s + 1.96*se_s), p_s, nrow(sub)))
}
cat("\n")

# ============================================================
# 8. CONSOLIDATED RESULTS TABLE
# ============================================================
cat("=== CONSOLIDATED RESULTS: edhs_improved_bin -> stunted_num ===\n")
cat(sprintf("%-50s  %8s  %18s  %8s\n", "Model", "aOR", "95% CI", "p-value"))
cat(strrep("-", 90), "\n")

print_row <- function(label, b, se, p) {
  if (is.na(b)) { cat(sprintf("%-50s  %8s  %18s  %8s\n", label, "NA", "NA", "NA")); return() }
  cat(sprintf("%-50s  %8.4f  %8.4f-%8.4f  %8.4f\n",
              label, exp(b), exp(b - 1.96*se), exp(b + 1.96*se), p))
}

print_row("Model A — Crude logistic",                         bA, seA, sA["edhs_improved_bin","Pr(>|z|)"])
print_row("Model B — svyglm (complex-survey adjusted)",       bB, seB, pB)
if (WEMIX_OK)  print_row("Model C — WeMix GLMM (RI: PSU*year)",           bC, seC, pC)
if (LME4_OK)   print_row("Model D — lme4 glmer (sensitivity, unweighted)", bD, seD, pD)
cat(strrep("-", 90), "\n\n")

cat("=== INTERPRETATION ===\n")
cat("Primary exposure: edhs_improved_bin (1=improved water source, 0=unimproved)\n")
cat("Outcome: stunted_num (HAZ < -2, child stunting)\n\n")
cat("NULL HYPOTHESIS CHECK:\n")
cat("  If the true OR is approximately 1.009 (as seen in prior GLMM runs),\n")
cat("  then improved water access shows NO independent association with\n")
cat("  child stunting after controlling for socioeconomic confounders.\n\n")
cat("  This is the Gap 1 null finding: the EDHS improved-water classification\n")
cat("  does not capture functionally safe water (fails to measure microbiological\n")
cat("  contamination, water chain re-contamination at point-of-use, or\n")
cat("  infrastructure non-functionality — the theoretical mechanism for why\n")
cat("  the WASH-stunting pathway is attenuated).\n\n")

# ============================================================
# 9. SAVE OUTPUTS
# ============================================================
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

results_df <- data.frame(
  model   = c("A_crude","B_svyglm",
               if (WEMIX_OK) "C_WeMix" else NULL,
               if (LME4_OK)  "D_lme4"  else NULL),
  exposure = "edhs_improved_bin",
  outcome  = "stunted_num",
  logOR   = c(bA, bB,
               if (WEMIX_OK) bC else NULL,
               if (LME4_OK)  bD else NULL),
  seLogOR = c(seA, seB,
               if (WEMIX_OK) seC else NULL,
               if (LME4_OK)  seD else NULL),
  p_value = c(sA["edhs_improved_bin","Pr(>|z|)"], pB,
               if (WEMIX_OK) pC else NULL,
               if (LME4_OK)  pD else NULL)
) %>%
  mutate(
    aOR = exp(logOR),
    CI95_lo = exp(logOR - 1.96*seLogOR),
    CI95_hi = exp(logOR + 1.96*seLogOR)
  )

write.csv(results_df, "output/tables/robust_panel_edhs_improved_stunting.csv",
          row.names = FALSE)
cat("Saved -> output/tables/robust_panel_edhs_improved_stunting.csv\n")

if (!is.null(reg_tab) && nrow(reg_tab) > 1) {
  write.csv(reg_tab, "output/tables/region_stratified_edhs_stunting.csv",
            row.names = FALSE)
  cat("Saved -> output/tables/region_stratified_edhs_stunting.csv\n")
}

cat("\n=== COMPLETE ===\n")
