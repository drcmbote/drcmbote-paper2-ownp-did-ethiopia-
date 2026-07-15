# Secular Trends in Child Stunting and Policy Attribution Under Ethiopia's One WASH National Programme

Reproducibility package for the analysis reported in:

> Secular Trends in Child Stunting and Policy Attribution Under Ethiopia's One WASH National Programme: A Six-Round Generalised Intensity Difference-in-Differences Analysis, 2000-2024. Submitted to *The Lancet Global Health*.

This repository contains the R and Python scripts used to construct the pooled six-round analytic panel and to produce every table and figure in the manuscript and its Supplementary Appendix. It does not contain raw data (see **Data access** below). This is Paper 2 of a two-paper series; Paper 1 (zone-level water point functionality and child nutritional outcomes, 2024-25 EDHS) has its own reproducibility repository.

## Data access

- **Ethiopian Demographic and Health Survey (EDHS)** individual child, woman, and household recode files for all six rounds used (2000, 2005, 2011, 2016, 2019, 2024-25) are available on registration at [dhsprogram.com](https://dhsprogram.com), subject to the DHS Program's standard data use agreement. Redistribution of the raw DHS extracts is not permitted under that agreement, which is why they are not included here.

Analysts who have obtained all six rounds under the DHS data use agreement can reproduce the pooled analytic panel (`panel_5rounds_child.csv`, referenced by the scripts below; retains its original five-round filename for continuity even though it was extended to six rounds) using `build_panel_5rounds.py`, which documents every round-specific data quirk handled (wealth-index and anthropometry merges for 2000/2005, `b19` reconstruction, region-code harmonisation across rounds, JMP improved/unimproved keyword classification, SNNPR post-2020 subdivision reunification, and the 2024-25 Gambella spelling fix).

## Script-to-output mapping

| Script | Produces |
|---|---|
| `build_panel_5rounds.py` | Pooled six-round analytic panel (`data/data_clean/panel_5rounds_child.csv`), the shared input for every script below |
| `table1_gap2.R` | Table 1 (baseline sociodemographic and WASH characteristics by survey round) |
| `crude_or_gap2.R` | Unadjusted/crude odds ratios (survey-round trend and main-effect crude OR reported in Results); feeds `figure1_p2.R` |
| `did_gap2_intensity.R` | Primary generalised intensity DiD model: Table 2 (survey-round fixed effects and OWNP policy interaction terms); feeds `figure1_p2.R` and `figure2_p2.R` |
| `did_gap2_sensitivity.R` | Model A (five-round, diarrhoea-adjusted, excludes 2019) and Model B (2000-baseline intensity): Table 3 Panel B |
| `event_study_placebo_gap2.R` | Event-study specification and placebo falsification test: Table 3 Panel A; feeds `figure2_p2.R` |
| `figure1_p2.R` | Figure 1 (observed stunting prevalence by round; adjusted vs crude survey-year odds ratios) |
| `figure2_p2.R` | Figure 2 (event-study coefficients; forest plot of the five DiD interaction terms) |
| `flowdiagram_p2.R` | Figure S1 (participant flow diagram, Supplementary Appendix) |
| `run_robust_panel_model.R` | Independent pooled main-effect robustness benchmark: crude/svyglm/lme4 models, region-stratified estimates, and the pre-/post-OWNP subgroup comparison reported in Table S1 (Supplementary Appendix) |

## Statistical environment

Models were fitted in R 4.5.1. The primary generalised intensity difference-in-differences model and its sensitivity/event-study/placebo variants used the `survey` package (`svyglm`, quasi-binomial link) with round-qualified PSU (`combined_psu = survey_year_v001`) and stratum (`combined_strata = survey_year_v023`) identifiers to prevent cross-round conflation of cluster/stratum codes, and DHS child sampling weights applied at the individual level. The independent main-effect robustness benchmark additionally used `lme4::glmer` (random intercept for `combined_psu`) and an unadjusted crude logistic regression, pooled across regions by inverse-variance fixed-effect meta-analysis (Cochran's Q, I²). Variance inflation factors were screened with the `car` package; all values were below 3.5. See the manuscript Methods for full specification details.

## License

Code is released under the MIT License. See `LICENSE`.
