#!/usr/bin/env python3
"""MIQE standard-curve and control QC, per plate x target.

The four control classes are gated separately because each implicates a
different failure and a different corrective action (assay design §5, §9):

    ntc               no amplification   PCR reagent contamination
    extraction_blank  no amplification   extraction reagent / consumable
                                         contamination
    negative_host     no amplification   cross-reactivity with the uninfected
                                         Anopheles genome
    positive          amplification      assay failure or inhibition, and
                                         separately a Cq outside tolerance of
                                         the run-wide per-target median

Collapsing the first three into one "negative control" count, as a single
gate does, reports a contaminated extraction and a cross-reacting primer set
as the same number.

Standard-curve acceptance follows assay design §4: R2 >= 0.99, efficiency
90-110%, slope -3.58 to -3.10.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from well_roles import (AMPLIFICATION_EXPECTED, CONTROL_ROLES,
                        NO_AMPLIFICATION_EXPECTED, POSITIVE)

ap = argparse.ArgumentParser()
ap.add_argument("--wells", required=True, type=Path)
ap.add_argument("--r2-min", type=float, default=0.99)
ap.add_argument("--eff-min", type=float, default=90.0)
ap.add_argument("--eff-max", type=float, default=110.0)
ap.add_argument("--slope-min", type=float, default=-3.58)
ap.add_argument("--slope-max", type=float, default=-3.10)
ap.add_argument("--ctl-amplification-max", type=float, default=40.0,
                help="a control well with a Cq at or below this has "
                     "amplified. Distinct from the calling rule's "
                     "--cq-trace-max, which bounds the trace band for "
                     "specimen wells; the two happen to share a value but "
                     "answer different questions")
ap.add_argument("--pos-ctl-cq-tolerance", type=float, default=2.0,
                help="cycles either side of the run-wide per-target median "
                     "positive-control Cq (assay design §9 expects "
                     "amplification within +/- 2 Ct of the historical mean; "
                     "absent a historical series, the run median is the "
                     "reference and is reported as such)")
ap.add_argument("--out-target", required=True, type=Path)
ap.add_argument("--out-control", required=True, type=Path)
a = ap.parse_args()

w = pd.read_csv(a.wells)

# ---------------------------------------------------------------- controls --
ctl = w[w.role.isin(CONTROL_ROLES)].copy()
amplified = ctl.cq.notna() & (ctl.cq <= a.ctl_amplification_max)
ctl["amplified"] = amplified
ctl["verdict"] = "pass"

for role, verdict in NO_AMPLIFICATION_EXPECTED.items():
    ctl.loc[ctl.role.eq(role) & amplified, "verdict"] = verdict
for role, verdict in AMPLIFICATION_EXPECTED.items():
    ctl.loc[ctl.role.eq(role) & ~amplified, "verdict"] = verdict

# Positive controls that amplified but drifted from the run-wide per-target
# median. Reported separately from outright failure: a drifting positive is a
# trend worth seeing before it becomes a failure.
pos = ctl.role.eq(POSITIVE) & amplified
if pos.any():
    ref = ctl.loc[pos].groupby("target").cq.median().rename("pos_ctl_cq_ref")
    ctl = ctl.merge(ref, left_on="target", right_index=True, how="left")
    drift = pos & (ctl.cq - ctl.pos_ctl_cq_ref).abs().gt(a.pos_ctl_cq_tolerance)
    ctl.loc[drift & ctl.verdict.eq("pass"), "verdict"] = "WARN_pos_ctl_cq_drift"
else:
    ctl["pos_ctl_cq_ref"] = pd.NA

ctl_out = ctl[["plate", "plate_class", "well", "sample_name", "role", "target",
               "control_target", "cq", "amplified", "pos_ctl_cq_ref",
               "verdict"]]
ctl_out.to_csv(a.out_control, index=False)

# Wells whose dedicated target disagrees with the well's Target Name would
# mean a control was pipetted into the wrong block.
mis = ctl[ctl.control_target.notna() & ctl.control_target.ne(ctl.target)]

# ------------------------------------------------------- standard curves ----
std = w[w.role == "standard"]
curve = (std.groupby(["plate", "plate_class", "target"])
            .agg(n_standard_wells=("cq", "size"),
                 n_levels=("quantity", "nunique"),
                 r2=("r2", "median"), slope=("slope", "median"),
                 efficiency=("efficiency", "median"),
                 y_intercept=("y_intercept", "median"))
            .reset_index())

# One failure count per control class, so the reason a plate x target failed
# is legible without going back to control_qc.csv.
fail = ctl.assign(bad=ctl.verdict.str.startswith("FAIL"))
counts = (fail.pivot_table(index=["plate", "target"], columns="role",
                           values="bad", aggfunc="sum", fill_value=0)
              .rename(columns={r: f"n_fail_{r}" for r in CONTROL_ROLES})
              .reset_index())
qc = curve.merge(counts, on=["plate", "target"], how="left")
for r in CONTROL_ROLES:
    col = f"n_fail_{r}"
    if col not in qc.columns:
        qc[col] = 0
    qc[col] = qc[col].fillna(0).astype(int)
qc["n_control_failures"] = qc[[f"n_fail_{r}" for r in CONTROL_ROLES]].sum(axis=1)

reasons = []
for _, r in qc.iterrows():
    bad = []
    if pd.isna(r.efficiency) or not (a.eff_min <= r.efficiency <= a.eff_max):
        bad.append(f"efficiency={r.efficiency}")
    if pd.isna(r.r2) or r.r2 < a.r2_min:
        bad.append(f"r2={r.r2}")
    if pd.notna(r.slope) and not (a.slope_min <= r.slope <= a.slope_max):
        bad.append(f"slope={r.slope}")
    for role in CONTROL_ROLES:
        n = int(r[f"n_fail_{role}"])
        if n:
            bad.append(f"{role}={n}")
    reasons.append(";".join(bad))
qc["fail_reasons"] = reasons
qc["curve_ok"] = qc.fail_reasons == ""
qc.to_csv(a.out_target, index=False)

print(f"plate x target curves: {int(qc.curve_ok.sum())}/{len(qc)} pass")
for role in CONTROL_ROLES:
    n = int(qc[f"n_fail_{role}"].sum())
    print(f"  {role:<17} failures: {n}")
if len(mis):
    print(f"  {len(mis)} control wells sit in a block for a different target:")
    print(mis[["plate", "well", "sample_name", "control_target", "target"]]
          .to_string(index=False))
