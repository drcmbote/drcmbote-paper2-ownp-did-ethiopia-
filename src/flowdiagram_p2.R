# Paper 2 (Gap 2) — Participant Flow Diagram
# STROBE Item 13(c): flow diagram showing pooled panel sample derivation
#   per survey round (2000, 2005, 2011, 2016, 2019)
# Outputs:
#   output/figures/flowdiagram_p2.tiff  (300 DPI, Lancet submission quality)
#   output/figures/flowdiagram_p2.png   (150 DPI preview)
#   output/tables/flowdiagram_p2_counts.csv

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
})

# =============================================================================
# STAGE 1: Count participants at each exclusion step
# =============================================================================
panel <- read.csv(
  "data/data_clean/panel_5rounds_child.csv",
  stringsAsFactors = FALSE, na.strings = c("NA", "")
)

rounds <- sort(unique(panel$survey_year))

# Stage A: total records across all 5 rounds
n_A_total <- nrow(panel)
n_A_round <- sapply(rounds, function(y) sum(panel$survey_year == y))

# Stage B: exclude missing stunting status
panel_b <- panel %>% filter(!is.na(stunted_num))
n_excl_B <- n_A_total - nrow(panel_b)

# Stage C: exclude missing child weight
panel_c <- panel_b %>% filter(!is.na(child_wt), child_wt > 0)
n_excl_C <- nrow(panel_b) - nrow(panel_c)

# Stage D: exclude missing improved water classification
panel_d <- panel_c %>% filter(!is.na(edhs_improved_bin))
n_excl_D <- nrow(panel_c) - nrow(panel_d)

# Stage E: intensity join (regions without 2011 baseline match)
intensity_tab <- panel %>%
  filter(survey_year == 2011, !is.na(edhs_improved_bin),
         !is.na(child_wt), child_wt > 0) %>%
  group_by(region) %>%
  summarise(intensity = sum((1 - edhs_improved_bin) * child_wt) / sum(child_wt),
            .groups = "drop")
panel_e <- panel_d %>% left_join(intensity_tab %>% select(region, intensity), by = "region") %>%
  filter(!is.na(intensity))
n_excl_E <- nrow(panel_d) - nrow(panel_e)

# Final per-round breakdown
n_final       <- nrow(panel_e)
n_final_round <- sapply(rounds, function(y) sum(panel_e$survey_year == y))
n_psu         <- length(unique(paste0(panel_e$survey_year, "_", panel_e$v001)))
n_strata      <- length(unique(paste0(panel_e$survey_year, "_", panel_e$v023)))

# Note about 2019: diarrhoea excluded from pooled model
n_2019_diarrhea_missing <- sum(panel_e$survey_year == 2019 &
                                is.na(panel_e$recent_diarrhea_num), na.rm = TRUE)

cat("=== FLOW DIAGRAM COUNTS (Paper 2) ===\n")
cat(sprintf("A  Total (%d-round EDHS, all regions):            N = %s\n",
            length(rounds), format(n_A_total, big.mark=",")))
for (i in seq_along(rounds))
  cat(sprintf("   Round %d:  n = %s\n", rounds[i], format(n_A_round[i], big.mark=",")))
cat(sprintf("   Excl. missing stunting status:               n = %s\n",
            format(n_excl_B, big.mark=",")))
cat(sprintf("   Excl. missing/zero sampling weight:          n = %s\n",
            format(n_excl_C, big.mark=",")))
cat(sprintf("   Excl. missing water source classification:   n = %s\n",
            format(n_excl_D, big.mark=",")))
cat(sprintf("   Excl. regions without 2011 intensity data:   n = %s\n",
            format(n_excl_E, big.mark=",")))
cat(sprintf("Final analytical sample:                         N = %s\n",
            format(n_final, big.mark=",")))
cat(sprintf("   %d round-qualified PSUs, %d strata\n", n_psu, n_strata))
for (i in seq_along(rounds))
  cat(sprintf("   Round %d: n = %s\n", rounds[i], format(n_final_round[i], big.mark=",")))
cat(sprintf("Note: 2019 interim DHS — diarrhoea (h11) absent in n = %s records\n",
            format(n_2019_diarrhea_missing, big.mark=",")))

counts <- data.frame(
  stage = c("A_total", paste0("A_round_", rounds),
            "excl_stunting", "excl_weight", "excl_water", "excl_intensity",
            "final_total", paste0("final_round_", rounds)),
  n     = c(n_A_total, n_A_round,
            n_excl_B, n_excl_C, n_excl_D, n_excl_E,
            n_final, n_final_round)
)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(counts, "output/tables/flowdiagram_p2_counts.csv", row.names = FALSE)

# =============================================================================
# STAGE 2: Build flow diagram with ggplot2
# =============================================================================
fmt <- function(...) paste0(sprintf(...))

# Layout: x in [0,10], y in [0,10]
# Main boxes centred at x=5; exclusion boxes right side; round boxes at bottom

n_rounds <- length(rounds)
gap <- 0.25
round_box_width <- (10 - gap * (n_rounds - 1)) / n_rounds

round_extra <- function(y) {
  if (y == 2011) return("\n(pre-OWNP\nreference)")
  if (y == 2019) return("\n(interim\nDHS)")
  ""
}
round_fill <- function(y) {
  if (y == 2011) return("#E8F4E8")
  if (y == 2019) return("#FFF8F0")
  "#F7F7F7"
}
round_fface <- function(y) if (y == 2011) "bold" else "plain"

round_boxes <- list()
for (i in seq_along(rounds)) {
  xl <- (i - 1) * (round_box_width + gap)
  round_boxes[[paste0("R", i)]] <- list(
    xl = xl, xr = xl + round_box_width, yb = 1.2, yt = 2.8,
    label = fmt("%d\nn = %s%s", rounds[i],
                format(n_final_round[i], big.mark = ","), round_extra(rounds[i]))
  )
}
round_cx <- sapply(round_boxes, function(b) (b$xl + b$xr) / 2)

boxes <- c(list(
  A = list(xl=1.5, xr=8.5, yb=8.8, yt=10,
           label = fmt("%d-round EDHS pooled dataset\n%s\nN = %s",
                       n_rounds, paste(rounds, collapse = " / "),
                       format(n_A_total, big.mark=","))),

  B = list(xl=1.5, xr=8.5, yb=5.5, yt=6.9,
           label = fmt("Analytical sample\n(complete stunting status, water classification,\nsampling weight, and regional intensity)\nN = %s — %d round-qualified PSUs, %d strata",
                       format(n_final, big.mark=","), n_psu, n_strata)),

  note = list(xl=1.5, xr=8.5, yb=4.2, yt=5.1,
              label = "Note: recent diarrhoea (h11) unavailable in 2019 interim DHS;\nexcluded from pooled model; re-included in Model A (sensitivity, all rounds except 2019)")
),
round_boxes,
list(
  # Exclusion box
  eAB = list(xl=7.5, xr=10,  yb=6.9, yt=8.5,
             label = fmt("Excluded\nMissing stunting status: n = %s\nMissing weight:          n = %s\nMissing water class.:    n = %s\nMissing region intensity:n = %s",
                         format(n_excl_B, big.mark=","),
                         format(n_excl_C, big.mark=","),
                         format(n_excl_D, big.mark=","),
                         format(n_excl_E, big.mark=",")))
))

draw_box <- function(p, b, fill="white", col="black", lwd=0.6) {
  p + annotate("rect",
               xmin=b$xl, xmax=b$xr, ymin=b$yb, ymax=b$yt,
               fill=fill, colour=col, linewidth=lwd)
}
add_label <- function(p, b, size=2.85, fontface="plain") {
  p + annotate("text",
               x=(b$xl+b$xr)/2, y=(b$yb+b$yt)/2,
               label=b$label, size=size, hjust=0.5, vjust=0.5,
               fontface=fontface, lineheight=1.15)
}

p <- ggplot() +
  scale_x_continuous(limits=c(0,10), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,10), expand=c(0,0)) +
  theme_void() +
  theme(plot.margin=margin(10,10,10,10))

fills <- c(list(A="white", B="white", note="#F0F4FF"),
           setNames(lapply(rounds, round_fill), paste0("R", seq_along(rounds))),
           list(eAB="#F7F7F7"))
fface <- c(list(A="plain", B="bold", note="italic"),
           setNames(lapply(rounds, round_fface), paste0("R", seq_along(rounds))),
           list(eAB="plain"))

for (nm in names(boxes)) {
  p <- draw_box(p, boxes[[nm]], fill=fills[[nm]])
  p <- add_label(p, boxes[[nm]], fontface=fface[[nm]])
}

arr <- arrow(length=unit(0.18,"cm"), type="closed")
mid_x <- 5

# A → B (vertical)
p <- p +
  annotate("segment", x=mid_x, xend=mid_x, y=8.8, yend=6.9,
           arrow=arr, linewidth=0.5)

# A → exclusion (horizontal tap-off)
p <- p +
  annotate("segment", x=mid_x, xend=mid_x, y=7.85, yend=7.85, linewidth=0) +
  annotate("segment", x=mid_x, xend=8.7, y=7.85, yend=7.85,
           arrow=arr, linewidth=0.4, linetype="dashed")

# B → note
p <- p +
  annotate("segment", x=mid_x, xend=mid_x, y=5.5, yend=5.1,
           arrow=arr, linewidth=0.4)

# note → split to rounds
split_y <- 3.9
p <- p +
  annotate("segment", x=mid_x, xend=mid_x, y=4.2, yend=split_y, linewidth=0.5)

for (cx in round_cx) {
  p <- p +
    annotate("segment", x=cx, xend=cx, y=split_y, yend=2.8,
             arrow=arr, linewidth=0.4)
}
p <- p +
  annotate("segment", x=round_cx[1], xend=round_cx[n_rounds], y=split_y, yend=split_y,
           linewidth=0.5)

# OWNP divider line between the last pre-OWNP round and the first post-OWNP round
is_post <- rounds >= 2016
last_pre_idx <- max(which(!is_post))
first_post_idx <- min(which(is_post))
divider_x <- (round_boxes[[last_pre_idx]]$xr + round_boxes[[first_post_idx]]$xl) / 2
p <- p +
  annotate("segment", x=divider_x, xend=divider_x, y=0.8, yend=3.2,
           linewidth=0.6, colour="grey40", linetype="dashed") +
  annotate("text", x=divider_x, y=0.55, label="← Pre-OWNP  |  Post-OWNP →",
           size=2.5, hjust=0.5, colour="grey40", fontface="italic")

# Title
p <- p +
  annotate("text", x=5, y=9.7,
           label="Figure. Participant flow diagram — Paper 2",
           size=3.4, fontface="bold", hjust=0.5)

# =============================================================================
# STAGE 3: Save
# =============================================================================
dir.create("output/figures", recursive=TRUE, showWarnings=FALSE)

ggsave("output/figures/flowdiagram_p2.tiff", p,
       width=8, height=9, dpi=300, compression="lzw")
cat("Saved -> output/figures/flowdiagram_p2.tiff\n")

ggsave("output/figures/flowdiagram_p2.png", p,
       width=8, height=9, dpi=150)
cat("Saved -> output/figures/flowdiagram_p2.png\n")
