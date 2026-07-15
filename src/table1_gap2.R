# Paper 2 (Gap 2) — Table 1: Baseline Characteristics by Survey Round
# STROBE Items 14(a) and 14(b): participant characteristics + missing data counts
# Columns: Variable | Category | one column per survey round | Missing n (%)
# Missing counts are pooled across all rounds (total panel missing)
# Output: output/tables/table1_gap2.csv

suppressMessages({
  library(dplyr)
  library(tidyr)
})

# ---- Load and prepare data (identical pipeline to did_gap2_intensity.R) ----
panel <- read.csv(
  "data/data_clean/panel_5rounds_child.csv",
  stringsAsFactors = FALSE, na.strings = c("NA", "")
)

intensity_tab <- panel %>%
  filter(survey_year == 2011, !is.na(edhs_improved_bin),
         !is.na(child_wt), child_wt > 0) %>%
  group_by(region) %>%
  summarise(intensity = sum((1 - edhs_improved_bin) * child_wt) / sum(child_wt),
            .groups = "drop")

df <- panel %>%
  left_join(intensity_tab %>% select(region, intensity), by = "region") %>%
  mutate(
    survey_year_f     = factor(survey_year),
    wealth_q          = factor(wealth_q,
                               levels = c("poorest","poorer","middle","richer","richest")),
    mother_education_f = factor(mother_education,
                                levels = c("no_education","primary","secondary","higher"),
                                labels = c("No education","Primary","Secondary","Higher"))
  ) %>%
  filter(!is.na(stunted_num), !is.na(child_wt), child_wt > 0,
         !is.na(edhs_improved_bin), !is.na(intensity))

rounds <- sort(unique(df$survey_year))
N_total <- nrow(df)
N_round <- sapply(rounds, function(y) sum(df$survey_year == y))

cat(sprintf("Table 1 analytical sample: N = %d\n", N_total))
cat("Per round:", paste(rounds, N_round, sep="=", collapse=" | "), "\n")

# ---- Helper functions -------------------------------------------------------
# DHS convention: n = unweighted case count, % = DHS-weighted proportion.
# (Previously this table reported an unweighted % despite the footnote
# claiming "n (weighted %)"; fixed 2026-07-15 after cross-checking against
# the DHS-weighted prevalence actually cited in Results/Abstract/Discussion.)
wpct <- function(indicator, weight) {
  valid <- !is.na(indicator) & !is.na(weight)
  n <- sum(indicator[valid], na.rm = TRUE)
  tot_w <- sum(weight[valid])
  if (tot_w == 0) return("—")
  w_pct <- 100 * sum(weight[valid & indicator]) / tot_w
  sprintf("%s (%.1f%%)", format(n, big.mark=","), w_pct)
}
mean_sd <- function(x, weight) {
  valid <- !is.na(x) & !is.na(weight)
  x <- x[valid]; w <- weight[valid]
  if (length(x) == 0) return("—")
  wmean <- sum(x * w) / sum(w)
  wsd <- sqrt(sum(w * (x - wmean)^2) / sum(w))
  sprintf("%.1f (%.1f)", wmean, wsd)
}
miss_cell <- function(vec, N) {
  n_miss <- sum(is.na(vec))
  if (n_miss == 0) "0 (0.0%)" else sprintf("%d (%.1f%%)", n_miss, 100*n_miss/N)
}

rd <- function(y) df %>% filter(survey_year == y)   # subset by round
Nr <- function(y) sum(df$survey_year == y)           # N for round y

make_row <- function(variable, category, vals_by_round, missing_pooled) {
  row <- data.frame(
    Variable = variable,
    Category = category,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(rounds))
    row[[paste0("Y", rounds[i])]] <- vals_by_round[i]
  row[["Missing_pooled"]] <- missing_pooled
  row
}

rows <- list()

# ---- 0. SAMPLE SIZE --------------------------------------------------------
rows[[1]] <- make_row(
  "Analytical sample (N)", "",
  sapply(rounds, function(y) format(Nr(y), big.mark=",")),
  "—"
)

# ---- 1. OUTCOME ------------------------------------------------------------
rows[[2]] <- make_row(
  "Child stunting (HAZ < −2)", "Yes",
  sapply(rounds, function(y) {
    r <- rd(y); wpct(r$stunted_num==1, r$child_wt)
  }),
  miss_cell(df$stunted_num, N_total)
)

# ---- 2. CHILD CHARACTERISTICS ----------------------------------------------
rows[[3]] <- make_row(
  "Age, months — mean (SD)", "",
  sapply(rounds, function(y) { r <- rd(y); mean_sd(r$child_age_months, r$child_wt) }),
  miss_cell(df$child_age_months, N_total)
)

rows[[4]] <- make_row(
  "Male sex", "Yes",
  sapply(rounds, function(y) {
    r <- rd(y); wpct(r$child_sex_male_num==1, r$child_wt)
  }),
  miss_cell(df$child_sex_male_num, N_total)
)

# ---- 3. WATER / WASH -------------------------------------------------------
rows[[5]] <- make_row(
  "Improved water source use", "Yes",
  sapply(rounds, function(y) {
    r <- rd(y); wpct(r$edhs_improved_bin==1, r$child_wt)
  }),
  miss_cell(df$edhs_improved_bin, N_total)
)

# ---- 4. HOUSEHOLD CHARACTERISTICS ------------------------------------------
wealth_lvls <- c("poorest","poorer","middle","richer","richest")
wealth_lbl  <- c("Poorest","Poorer","Middle","Richer","Richest")
for (i in seq_along(wealth_lvls)) {
  lvl <- wealth_lvls[i]
  rows[[length(rows)+1]] <- make_row(
    if (i==1) "Wealth quintile" else "",
    wealth_lbl[i],
    sapply(rounds, function(y) {
      r <- rd(y); wpct(r$wealth_q==lvl, r$child_wt)
    }),
    if (i==1) miss_cell(df$wealth_q, N_total) else ""
  )
}

# ---- 5. MATERNAL EDUCATION -------------------------------------------------
edu_lvls <- c("No education","Primary","Secondary","Higher")
for (i in seq_along(edu_lvls)) {
  lvl <- edu_lvls[i]
  rows[[length(rows)+1]] <- make_row(
    if (i==1) "Maternal education" else "",
    lvl,
    sapply(rounds, function(y) {
      r <- rd(y); wpct(r$mother_education_f==lvl, r$child_wt)
    }),
    if (i==1) miss_cell(df$mother_education_f, N_total) else ""
  )
}

# ---- 6. OWNP PERIOD NOTE ---------------------------------------------------
intensity_baseline_year <- 2011
rows[[length(rows)+1]] <- make_row(
  "OWNP period", "",
  sapply(rounds, function(y) {
    post <- unique(df$post_ownp[df$survey_year == y])
    if (y == intensity_baseline_year) "Pre (reference)"
    else if (post == 0) "Pre" else "Post"
  }),
  "—"
)

# ---- 7. NOTE: diarrhoea ----------------------------------------------------
rows[[length(rows)+1]] <- make_row(
  "Recent diarrhoea (h11)", "Available",
  sapply(rounds, function(y) if (all(is.na(rd(y)$recent_diarrhea_num))) "Not fielded" else "Yes"),
  sprintf("%d (%.1f%%)", sum(is.na(df$recent_diarrhea_num)), 100*mean(is.na(df$recent_diarrhea_num)))
)

# ---- Compile ---------------------------------------------------------------
table1 <- do.call(rbind, rows)

yr_hdrs <- sprintf("Round %d\n(n = %s)", rounds, format(N_round, big.mark=","))
header <- data.frame(
  Variable = sprintf("Table 1. Baseline characteristics by survey round (N = %s)", format(N_total, big.mark=",")),
  Category = "",
  stringsAsFactors = FALSE
)
for (i in seq_along(rounds)) header[[paste0("Y",rounds[i])]] <- yr_hdrs[i]
header[["Missing_pooled"]] <- "Missing (pooled), n (%)"

table1_out <- rbind(header, table1)

write.csv(table1_out, "output/tables/table1_gap2.csv", row.names=FALSE, na="")
cat("Saved -> output/tables/table1_gap2.csv\n")
print(table1, row.names=FALSE)
