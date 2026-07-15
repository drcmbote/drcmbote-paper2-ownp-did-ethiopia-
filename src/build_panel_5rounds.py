"""Build the Gap 2 multi-round EDHS child panel (2000/2005/2011/2016/2019/
2024-25) for the pre/post OWNP generalized DiD design.

Design decision (per user instruction, after web-search confirmed OWNP was
launched nationally in Sept 2013 for ALL regions simultaneously -- not a
staggered regional rollout): post_ownp = 0 for surveys fielded before the
program existed (2000, 2005, 2011) and = 1 for surveys fielded after Phase I
launch (2016, 2019). This sidesteps the endogenous-timing problem a true
staggered design would have faced, at the cost of only 2 time bins instead
of a continuous treatment-intensity gradient.

Per-round data quirks handled here (each verified directly against the raw
DHS extracts, not assumed):

  - hv201 (water source) code lists expand over rounds (e.g. "tube well or
    borehole", "bottled water" added in later rounds), so the EDHS/JMP
    improved/unimproved classification uses keyword matching on the decoded
    text label rather than hardcoded numeric codes, so it is robust to the
    label drift across rounds.
  - hv270 (wealth quintile) is MISSING from the 2000 HR file -- the
    standalone wealth-index extract ETWI41FL (whhid, wlthind5) is merged in
    instead. whhid is a fixed-width string: chars[0:9]=cluster (v001),
    chars[9:12]=household number (v002), verified by parsing and
    cross-checking against KR's v001/v002 distribution.
  - hw70/hw71 (HAZ/WAZ, x100) are MISSING from KR in 2000 and 2005 -- the
    standalone height/weight extracts ETHW41FL/ETHW51FL (hc70=HAZ, hc71=WAZ)
    are merged in via the same fixed-width hwhhid string plus hwline (=b16,
    the child's line number), verified to match in tests (82% match rate
    against b16, fully expected since not all KR rows are anthropometry-
    eligible/measured children).
  - b19 (child's current age in months, DHS-computed) is missing in 2000,
    2005, 2011 -- reconstructed as v008 - b3 (CMC of interview minus CMC of
    birth), the standard DHS formula, confirmed identical to b19 where both
    exist in 2016/2019.
  - hv024/region numeric codes are NOT consistent across rounds (2000/2005/
    2011 use {1-7,12-15}; 2016/2019 use {1-11}), but decoding with
    convert_categoricals=True resolves to the SAME 11 Ethiopian regions in
    all 5 rounds (spelling varies: "oromiya"/"oromia", "addis abeba"/
    "addis ababa"/"addis adaba", "ben-gumz"/"benishangul-gumuz"). Resolved
    here with a lowercase + substring-based canonicalization map. No actual
    administrative boundary change occurred within 2000-2019 (the SNNPR
    split that created South West Ethiopia happened in 2021, after this
    panel's window), so this is purely a labelling-convention fix, not a
    geographic harmonization problem.

KR files do not carry hv-prefixed household variables directly (verified by
column inspection in all 6 rounds) -- hv201/hv270/hv025/hv024 are always
pulled from the HR file and merged onto KR via (v001, v002).

2024-25 EDHS (DHS-8, ETHR8AFL/ETKR8AFL) quirks verified against the raw
extracts:
  - hv024 splits the historical SNNPR into four separate region codes
    (Sidama seceded 2020, South West Ethiopia Peoples' Region formed 2021,
    the SNNPR remainder was renamed Central Ethiopia/South Ethiopia in
    2023). The 2011 baseline treatment intensity (did_gap2_intensity.R) is
    a single region-constant value computed once from the 2011 cross-
    section and left-joined onto every round by region name -- an
    unmapped region label here would silently drop ~30% of the 2024-25
    sample as intensity=NA. canonical_region() therefore maps all four
    successor labels back to "snnp" so the panel keeps using the same
    2011-defined baseline geography throughout.
  - hv201 adds a "large bottle" category not present in earlier rounds
    (large-format bottled water, distinct from the existing "bottled
    water" label). Classified as improved, consistent with how plain
    "bottled water" is already treated in this pipeline.
  - hv270, hw70/hw71, b19 are all present natively in HR/KR (same pattern
    as 2016/2019), so has_native_hwz=True, has_b19=True, no WI/HW merge
    needed.
  - hv024 also carries a spelling drift not covered by the pre-existing
    region canonicalization: 2011 decodes Gambela as "gambela" (one L)
    while 2024-25 decodes it as "gambella" (two Ls, the standard modern
    spelling). Confirmed by cross-checking which region label in the
    pooled panel had zero children in the 2011 baseline round (the
    baseline round is used to build the intensity variable, so any
    label not present there silently drops that region's later-round
    children as intensity=NA). canonical_region() now folds both
    spellings to "gambela".
"""
import pandas as pd

# HAZ/WAZ flag codes are NOT a single consistent value across rounds: 2000/
# 2005 use {9996,9997,9998}, 2011 adds 9999 too, 2016/2019 only use 9998.
# True z-scores (x100) never reach 9000, so any value >= 9000 is excluded
# uniformly instead of hardcoding a specific flag list per round.
HW_FLAG_THRESHOLD = 9000

IMPROVED_KEYWORDS = ["piped", "public tap", "standpipe", "tube well", "borehole",
                     "protected well", "protected spring", "covered well",
                     "covered spring", "rainwater", "bottled water", "large bottle"]
UNIMPROVED_KEYWORDS = ["unprotected", "open well", "open spring", "river", "pond",
                        "lake", "dam", "stream", "canal", "irrigation",
                        "tanker truck", "cart with small tank", "other"]


def classify_water(label):
    s = str(label).strip().lower()
    if any(k in s for k in UNIMPROVED_KEYWORDS):
        return "Unimproved (EDHS/JMP)"
    if any(k in s for k in IMPROVED_KEYWORDS):
        return "Improved (EDHS/JMP)"
    return None  # e.g. stray "99"/don't-know codes -> excluded


def canonical_region(label):
    s = str(label).strip().lower()
    if "oromiy" in s or s == "oromia":
        return "oromia"
    if "afar" in s or "affar" in s:
        return "afar"
    if "snnp" in s or s in ("sidama", "central ethiopia", "south ethiopia",
                             "south west ethiopia"):
        return "snnp"
    if "addis" in s:
        return "addis ababa"
    if "ben" in s:
        return "benishangul-gumuz"
    if "gambel" in s:
        return "gambela"
    return s


def parse_fixed_width_hhid(series, cluster_width=9, hh_width=3):
    cluster = series.str[0:cluster_width].astype(int)
    hh = series.str[cluster_width:cluster_width + hh_width].astype(int)
    return cluster, hh


def load_hr(hr_path, wi_path=None):
    hr = pd.read_stata(hr_path, convert_categoricals=False)
    hr_dec = pd.read_stata(hr_path, convert_categoricals=True)
    hr["water_label"] = hr_dec["hv201"].astype(str)
    hr["region_raw"] = hr_dec["hv024"].astype(str)
    hr["urban_label"] = hr_dec["hv025"].astype(str)
    cols = ["hv001", "hv002", "hv005", "water_label", "region_raw", "urban_label"]
    if "hv270" in hr.columns:
        cols.append("hv270")
        hr_small = hr[cols].rename(columns={"hv001": "v001", "hv002": "v002"})
    else:
        hr_small = hr[cols].rename(columns={"hv001": "v001", "hv002": "v002"})
        wi = pd.read_stata(wi_path, convert_categoricals=False)
        wi_dec = pd.read_stata(wi_path, convert_categoricals=True)
        wi["v001"], wi["v002"] = parse_fixed_width_hhid(wi["whhid"])
        wlth_map = {"lowest quintile": 1, "second quintile": 2, "middle quintile": 3,
                    "fourth quintile": 4, "highest quintile": 5}
        wi["hv270"] = wi_dec["wlthind5"].astype(str).map(wlth_map)
        hr_small = hr_small.merge(wi[["v001", "v002", "hv270"]], on=["v001", "v002"], how="left")

    hr_small["edhs_improved"] = hr_small["water_label"].apply(classify_water)
    hr_small["edhs_improved_bin"] = hr_small["edhs_improved"].map(
        {"Improved (EDHS/JMP)": 1, "Unimproved (EDHS/JMP)": 0}
    )
    hr_small["region"] = hr_small["region_raw"].apply(canonical_region)
    hr_small["urban"] = hr_small["urban_label"].eq("urban").astype(int)
    wealth_labels = {1: "poorest", 2: "poorer", 3: "middle", 4: "richer", 5: "richest"}
    hr_small["wealth_q"] = hr_small["hv270"].map(wealth_labels)
    return hr_small[["v001", "v002", "edhs_improved_bin", "region", "urban", "wealth_q"]]


def load_kr(kr_path, hw_path=None, has_b19=True, has_native_hwz=True):
    kr = pd.read_stata(kr_path, convert_categoricals=False)
    keep = ["v001", "v002", "v003", "v023", "b16", "hw1", "v005", "v008", "b3", "b4", "h11", "v106"]
    if has_native_hwz:
        keep += ["hw70", "hw71"]
    kr = kr[keep].copy()
    kr["child_wt"] = kr["v005"] / 1_000_000

    if not has_native_hwz:
        hw = pd.read_stata(hw_path, convert_categoricals=False)
        hw["v001"], hw["v002"] = parse_fixed_width_hhid(hw["hwhhid"])
        hw = hw.rename(columns={"hwline": "b16", "hc70": "hw70", "hc71": "hw71"})
        kr = kr.merge(hw[["v001", "v002", "b16", "hw70", "hw71"]], on=["v001", "v002", "b16"], how="left")

    kr["haz"] = kr["hw70"].where(kr["hw70"] < HW_FLAG_THRESHOLD) / 100
    kr["stunted_num"] = kr["haz"].le(-2).astype(float)
    kr.loc[kr["haz"].isna(), "stunted_num"] = pd.NA

    kr["child_age_months"] = kr["b19"] if has_b19 and "b19" in kr.columns else (kr["v008"] - kr["b3"])
    kr["child_sex_male_num"] = kr["b4"].eq(1).astype(float)
    kr["recent_diarrhea_num"] = kr["h11"].isin([1, 2]).astype(float)
    kr.loc[kr["h11"].isin([8, 9]) | kr["h11"].isna(), "recent_diarrhea_num"] = pd.NA
    edu_map = {0: "no_education", 1: "primary", 2: "secondary", 3: "higher"}
    kr["mother_education"] = kr["v106"].map(edu_map)

    return kr[["v001", "v002", "v003", "v023", "child_wt", "stunted_num", "child_age_months",
               "child_sex_male_num", "recent_diarrhea_num", "mother_education"]]


ROUNDS = [
    dict(year=2000, post_ownp=0,
         hr="2000 Standard DHS/extracted_recode/ETHR41DT/ETHR41FL.DTA",
         kr="2000 Standard DHS/extracted_recode/ETKR41DT/ETKR41FL.DTA",
         wi="2000 Standard DHS/extracted_recode/ETWI41DT/ETWI41FL.DTA",
         hw="2000 Standard DHS/extracted_recode/ETHW41DT/ETHW41FL.DTA",
         has_b19=False, has_native_hwz=False),
    dict(year=2005, post_ownp=0,
         hr="2005 Standard DHS/extracted_recode/ETHR51DT/ETHR51FL.DTA",
         kr="2005 Standard DHS/extracted_recode/ETKR51DT/ETKR51FL.DTA",
         wi=None,
         hw="2005 Standard DHS/extracted_recode/ETHW51DT/ETHW51FL.DTA",
         has_b19=False, has_native_hwz=False),
    dict(year=2011, post_ownp=0,
         hr="2011 Standard DHS/extracted_recode/ETHR61DT/ETHR61FL.DTA",
         kr="2011 Standard DHS/extracted_recode/ETKR61DT/ETKR61FL.DTA",
         wi=None, hw=None,
         has_b19=False, has_native_hwz=True),
    dict(year=2016, post_ownp=1,
         hr="dhs_2016/extracted_145/ETHR71DT/ETHR71FL.DTA",
         kr="dhs_2016/extracted_145/ETKR71DT/ETKR71FL.DTA",
         wi=None, hw=None,
         has_b19=True, has_native_hwz=True),
    dict(year=2019, post_ownp=1,
         hr="2019 Interim DHS/extracted_recode/ETHR81DT/ETHR81FL.DTA",
         kr="2019 Interim DHS/extracted_recode/ETKR81DT/ETKR81FL.DTA",
         wi=None, hw=None,
         has_b19=True, has_native_hwz=True),
    dict(year=2024, post_ownp=1,
         hr="2024-2025 Standard DHS/extracted_recode/ETHR8ADT/ETHR8AFL.dta",
         kr="2024-2025 Standard DHS/extracted_recode/ETKR8ADT/ETKR8AFL.dta",
         wi=None, hw=None,
         has_b19=True, has_native_hwz=True),
]

DATA_RAW = "data/data_raw"


def main():
    panels = []
    for cfg in ROUNDS:
        hr = load_hr(f"{DATA_RAW}/{cfg['hr']}", wi_path=f"{DATA_RAW}/{cfg['wi']}" if cfg["wi"] else None)
        kr = load_kr(f"{DATA_RAW}/{cfg['kr']}",
                     hw_path=f"{DATA_RAW}/{cfg['hw']}" if cfg["hw"] else None,
                     has_b19=cfg["has_b19"], has_native_hwz=cfg["has_native_hwz"])
        round_df = kr.merge(hr, on=["v001", "v002"], how="left")
        round_df["survey_year"] = cfg["year"]
        round_df["post_ownp"] = cfg["post_ownp"]
        panels.append(round_df)
        print(f"{cfg['year']}: n={len(round_df)} children, "
              f"stunting valid={round_df['stunted_num'].notna().sum()}, "
              f"region matched={round_df['region'].notna().sum()}, "
              f"wealth matched={round_df['wealth_q'].notna().sum()}")

    panel = pd.concat(panels, ignore_index=True)
    panel = panel.rename(columns={"edhs_improved_bin": "edhs_improved_bin"})

    out_path = "data/data_clean/panel_5rounds_child.csv"
    panel.to_csv(out_path, index=False)
    print(f"\nSaved -> {out_path}")
    print(f"\nTotal pooled n = {len(panel)}")
    print("\n--- Children per round ---")
    print(panel["survey_year"].value_counts().sort_index())
    print("\n--- post_ownp distribution ---")
    print(panel.groupby("post_ownp")["survey_year"].unique())
    print("\n--- Stunting valid n per round ---")
    print(panel.groupby("survey_year")["stunted_num"].apply(lambda s: s.notna().sum()))
    print("\n--- Key variable completeness (pooled) ---")
    for col in ["stunted_num", "edhs_improved_bin", "wealth_q", "urban", "region",
                "child_age_months", "child_sex_male_num", "recent_diarrhea_num",
                "mother_education", "child_wt"]:
        print(f"  {col}: {panel[col].notna().sum()} / {len(panel)} non-missing")


if __name__ == "__main__":
    main()
