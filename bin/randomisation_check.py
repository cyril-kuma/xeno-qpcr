#!/usr/bin/env python3
"""Did the randomisation work, and is detection still plate-structured?

Two questions, and keeping them apart is the point of this stage.

1. DECORRELATION. The randomisation (assay design §8, executed by
   scripts/randomise.py at seed 20261001) exists so that plate is not
   confounded with collection site, bioclimatic zone, collection method or
   year. Whether it succeeded is a property of the allocation alone and is
   checkable before a single plate is run: a chi-square test of independence
   between plate and each stratifying variable, plus the largest share any one
   level takes on any one plate. Large p-values and modest maximum shares mean
   the allocation is balanced.

2. RESIDUAL PLATE STRUCTURE. Given a verified allocation, detection that still
   clusters by plate cannot be zone confounding -- the allocation rules that
   out. It is a technical plate effect: a reagent lot, a run, an operator, a
   thermal block. The permutation test counts plates sitting at exactly 0% or
   100% positivity against a null that shuffles outcomes across plates of the
   same sizes.

   The earlier version of this check reported clustering and concluded that
   prevalence was uninterpretable. With the allocation verified, the right
   conclusion is narrower and stronger: the plate term belongs in the model as
   a random effect or a cluster variable, and zone contrasts estimated within
   plates remain interpretable.

   Intraclass correlation makes the comparison concrete. Reporting the ICC of
   detection by plate next to the ICC by site and by zone shows which level
   the variance actually sits at. A plate ICC far above the zone ICC is a
   technical signal; the reverse is biology.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

ap = argparse.ArgumentParser()
ap.add_argument("--mosquito", required=True, type=Path)
ap.add_argument("--wells", required=True, type=Path)
ap.add_argument("--assignment", required=True, type=Path)
ap.add_argument("--specimens", required=True, type=Path,
                help="the specimens tab; the assignment records position and "
                     "identity only, so the attributes the allocation was "
                     "balanced on are joined from here on sample_code")
ap.add_argument("--strata", default="bioclimatic_zone,collection_method,"
                                    "collection_location,year",
                help="specimen attributes tested for independence from plate")
ap.add_argument("--icc-groups", default="plate,collection_location,"
                                        "bioclimatic_zone",
                help="grouping variables whose ICC of detection is reported")
ap.add_argument("--max-level-share", type=float, default=0.75,
                help="a single level taking more than this share of any one "
                     "plate is flagged as an imbalance")
ap.add_argument("--n-perm", type=int, default=2000)
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--out-diagnostic", required=True, type=Path)
ap.add_argument("--out-plate", required=True, type=Path)
ap.add_argument("--out-balance", required=True, type=Path)
a = ap.parse_args()


def icc_oneway(values: np.ndarray, groups: np.ndarray) -> float:
    """One-way random-effects intraclass correlation, unbalanced groups.

    The ANOVA estimator: between- and within-group mean squares, with the
    harmonic-style group-size correction k0 that unbalanced designs need. A
    negative variance estimate is reported as 0, which is what it means --
    no detectable clustering at this level.
    """
    df = pd.DataFrame({"y": values, "g": groups}).dropna()
    if df.g.nunique() < 2 or len(df) <= df.g.nunique():
        return np.nan
    grand = df.y.mean()
    agg = df.groupby("g").y.agg(["size", "mean"])
    k, n = len(agg), len(df)
    ssb = float((agg["size"] * (agg["mean"] - grand) ** 2).sum())
    ssw = float(((df.y - df.g.map(agg["mean"])) ** 2).sum())
    msb, msw = ssb / (k - 1), ssw / (n - k)
    k0 = (n - (agg["size"] ** 2).sum() / n) / (k - 1)
    denom = msb + (k0 - 1) * msw
    if denom <= 0:
        return 0.0
    return float(max(0.0, (msb - msw) / denom))


m = pd.read_csv(a.mosquito)
w = pd.read_csv(a.wells)
asg = pd.read_csv(a.assignment)
spec = pd.read_csv(a.specimens)

if "plate_class" in w.columns:
    w = w[w.plate_class.eq("routine")]
pm = (w[w.role.eq("sample")][["sample_code", "plate"]]
      .drop_duplicates("sample_code"))
d = m.merge(pm, on="sample_code", how="inner")

# ------------------------------------------------- 1. decorrelation check --
# The allocation is tested as allocated, over every slot, including the
# specimens that have no record. Restricting to the analysed set first would
# test the balance of a subset rather than of the randomisation.
join_cols = [c for c in spec.columns
             if c == "sample_code" or c not in asg.columns]
asg = asg.merge(spec[join_cols], on="sample_code", how="left")
if "collection_date" in asg.columns and "year" not in asg.columns:
    asg["year"] = asg.collection_date.astype(str).str[-4:]
strata = [s.strip() for s in a.strata.split(",") if s.strip()]
bal_rows = []
for s in strata:
    if s not in asg.columns:
        bal_rows.append({"variable": s, "test": "chi2_independence",
                         "statistic": np.nan, "p_value": np.nan,
                         "max_level_share": np.nan, "verdict": "column absent"})
        continue
    t = pd.crosstab(asg.plate_id, asg[s].fillna("no metadata"))
    chi2, p, _dof, _exp = chi2_contingency(t)
    share = float(t.div(t.sum(axis=1), axis=0).max().max())
    bal_rows.append({
        "variable": s, "test": "chi2_independence",
        "statistic": round(float(chi2), 2), "p_value": round(float(p), 4),
        "max_level_share": round(share, 3),
        "verdict": ("balanced" if p >= 0.05 and share <= a.max_level_share
                    else "IMBALANCED")})
bal = pd.DataFrame(bal_rows)
bal.to_csv(a.out_balance, index=False)

decorrelated = bool((bal.verdict == "balanced").all())

# ------------------------------------- 2. residual plate structure check --
g = d.groupby("plate").any_plasmodium_detected.agg(["size", "sum"])
g["pct_positive"] = (100 * g["sum"] / g["size"]).round(1)
g.reset_index().to_csv(a.out_plate, index=False)

obs = int(((g.pct_positive == 0) | (g.pct_positive == 100)).sum())
rng = np.random.default_rng(a.seed)
vals = d.any_plasmodium_detected.to_numpy().astype(float)
sizes = g["size"].to_numpy()
null = np.empty(a.n_perm, dtype=int)
for i in range(a.n_perm):
    sh = rng.permutation(vals)
    idx = ext = 0
    for s in sizes:
        blk = sh[idx:idx + s]
        idx += s
        if blk.mean() in (0.0, 1.0):
            ext += 1
    null[i] = ext
p_perm = float((null >= obs).mean())
clustered = bool(p_perm < 0.001 and obs > 2 * max(np.median(null), 1))

iccs = {}
for col in [c.strip() for c in a.icc_groups.split(",") if c.strip()]:
    if col in d.columns:
        iccs[col] = icc_oneway(vals, d[col].to_numpy())

if not decorrelated:
    flag = "CONFOUNDED"
    note = ("the allocation itself is imbalanced, so plate and the flagged "
            "variable cannot be separated; resolve the allocation before "
            "interpreting any plate-level result")
elif clustered:
    flag = "PLATE_EFFECT"
    note = ("the allocation is balanced, so plate-clustered detection is a "
            "technical effect rather than zone confounding; carry plate as a "
            "random effect or cluster variable and estimate zone contrasts "
            "within plates")
else:
    flag = "ok"
    note = ("the allocation is balanced and detection shows no plate "
            "clustering beyond chance")

rows = [
    {"statistic": "randomisation_seed_recorded", "value": "20261001"},
    {"statistic": "n_plates", "value": int(len(g))},
    {"statistic": "n_specimens", "value": int(len(d))},
    {"statistic": "allocation_decorrelated", "value": decorrelated},
    {"statistic": "overall_positivity_pct", "value": round(100 * vals.mean(), 1)},
    {"statistic": "n_plates_all_negative", "value": int((g.pct_positive == 0).sum())},
    {"statistic": "n_plates_all_positive", "value": int((g.pct_positive == 100).sum())},
    {"statistic": "n_plates_extreme_observed", "value": obs},
    {"statistic": "n_plates_extreme_expected_median", "value": float(np.median(null))},
    {"statistic": "permutation_p_value", "value": p_perm},
]
for col, v in iccs.items():
    rows.append({"statistic": f"icc_detection_by_{col}",
                 "value": round(v, 4) if pd.notna(v) else np.nan})
rows += [{"statistic": "batch_effect_flag", "value": flag},
         {"statistic": "interpretation", "value": note}]

out = pd.DataFrame(rows)
out.to_csv(a.out_diagnostic, index=False)

print("\nallocation balance")
print(bal.to_string(index=False))
print("\nresidual plate structure")
print(out.to_string(index=False))
