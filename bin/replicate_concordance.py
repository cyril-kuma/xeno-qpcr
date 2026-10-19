#!/usr/bin/env python3
"""Inter-plate replicate concordance.

The randomisation left six slots spare and filled them with a second aliquot
of a specimen already placed on a different plate (the sample sheet's
`assignment` tab, role `INTER-PLATE REPLICATE (of plate Rnn)`). Those six
pairs are the only measurement in the design that separates plate-to-plate
technical variation from specimen-to-specimen biological variation: the same
DNA, the same assay, two plates.

Six pairs is a small number, and this stage reports them rather than testing
them. Per pair and per target it gives both calls, both Cq values and the
difference, plus whether the pair agrees on the call. The aggregate agreement
rate and the mean absolute Cq difference are the figures that belong next to
the batch diagnostic: a batch effect that moves calls should show up here as
disagreement, and one that does not is a weaker claim.

INDEPENDENCE
------------
Two reactions in two wells on two plates do not produce the same Cq to six
decimal places, and they do not produce the same Cq-confidence either. When
they do, the second record is a copy of the first rather than a second
measurement, and the pair says nothing about reproducibility however well it
agrees. Each pair is therefore checked for independence before its agreement
is counted, and a pair that fails is reported as `not_independent` and
excluded from the aggregate. Agreement computed over copies would read as
perfect concordance and would be the most misleading number in the report.

Specimens also carry `first_run_plate`, the plate they sat on in the first
round of this assay. Supplying that round's calls with `--round1-calls` adds a
second, larger concordance comparison across the whole set. Without it the
column is reported as present but uncompared, since the first-round result
table is not part of this repository.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--calls", required=True, type=Path,
                help="well-level target calls from call_targets.py")
ap.add_argument("--assignment", required=True, type=Path)
ap.add_argument("--targets", default="Pf,Pm,Po")
ap.add_argument("--round1-calls", type=Path,
                help="optional first-round mosquito-level table with "
                     "sample_code and <target>_call columns")
ap.add_argument("--independence-eps", type=float, default=1e-6,
                help="two wells whose Cq and Cq-confidence agree to within "
                     "this are treated as the same record duplicated rather "
                     "than as two measurements")
ap.add_argument("--out-pairs", required=True, type=Path)
ap.add_argument("--out-summary", required=True, type=Path)
a = ap.parse_args()

targets = [t.strip() for t in a.targets.split(",") if t.strip()]
calls = pd.read_csv(a.calls)
asg = pd.read_csv(a.assignment)

calls["plate"] = calls.plate.astype(str).str.upper()
calls["well"] = calls.well.astype(str).str.upper()
asg["plate_id"] = asg.plate_id.astype(str).str.upper()

reps = asg[asg.is_replicate.astype(str).str.lower().isin(["true", "1"])].copy()
rows = []

for _, r in reps.iterrows():
    code = str(r.sample_code).replace("_", "").upper()
    origin = (str(r.replicate_of_plate).upper()
              if pd.notna(r.replicate_of_plate) else None)
    for t in targets:
        well = str(r[f"well_{t}"]).upper()
        rep = calls[(calls.plate == r.plate_id) & (calls.well == well)
                    & (calls.target == t)]
        # The primary aliquot is the same specimen code on any other plate.
        prim = calls[(calls.sample_code.astype(str).str.upper() == code)
                     & (calls.target == t) & (calls.plate != r.plate_id)]
        if origin:
            on_origin = prim[prim.plate == origin]
            if len(on_origin):
                prim = on_origin
        if rep.empty or prim.empty:
            rows.append({"sample_code": r.sample_code, "target": t,
                         "replicate_plate": r.plate_id,
                         "primary_plate": origin, "status": "incomplete"})
            continue
        rep, prim = rep.iloc[0], prim.iloc[0]
        d = (rep.cq - prim.cq if pd.notna(rep.cq) and pd.notna(prim.cq)
             else np.nan)

        # Independence: the two records must differ somewhere. Identical Cq
        # and identical Cq-confidence means one is a copy of the other, and
        # the pair carries no information about plate-to-plate variation.
        same_cq = (pd.isna(rep.cq) and pd.isna(prim.cq)) or (
            pd.notna(d) and abs(d) <= a.independence_eps)
        conf_r, conf_p = rep.get("cq_conf"), prim.get("cq_conf")
        same_conf = (pd.isna(conf_r) and pd.isna(conf_p)) or (
            pd.notna(conf_r) and pd.notna(conf_p)
            and abs(float(conf_r) - float(conf_p)) <= a.independence_eps)
        # Two wells that both returned no amplification agree trivially and
        # carry no Cq to compare, so they are uninformative rather than
        # evidence that the record was duplicated.
        no_data = pd.isna(rep.cq) and pd.isna(prim.cq)
        independent = not (same_cq and same_conf) or no_data

        rows.append({
            "sample_code": r.sample_code, "target": t,
            "primary_plate": prim.plate, "primary_well": prim.well,
            "primary_call": prim.call, "primary_cq": prim.cq,
            "primary_cq_conf": conf_p,
            "replicate_plate": r.plate_id, "replicate_well": rep.well,
            "replicate_call": rep.call, "replicate_cq": rep.cq,
            "replicate_cq_conf": conf_r,
            "delta_cq": round(d, 3) if pd.notna(d) else np.nan,
            "call_agrees": bool(rep.call == prim.call),
            "independent": independent,
            "stratum": r.get("stratum"),
            "status": ("compared" if independent
                       else "not_independent" if not no_data
                       else "compared")})

pairs = pd.DataFrame(rows)
pairs.to_csv(a.out_pairs, index=False)

n_dep = int((pairs.status == "not_independent").sum()) if len(pairs) else 0
done = pairs[pairs.status == "compared"] if len(pairs) else pairs
n_cmp = len(done)
n_agree = int(done.call_agrees.sum()) if n_cmp else 0
mad = (float(done.delta_cq.abs().mean())
       if n_cmp and done.delta_cq.notna().any() else np.nan)
n_informative = (int(done.delta_cq.notna().sum()) if n_cmp else 0)

summary = [
    {"statistic": "n_replicate_specimens",
     "value": int(reps.sample_code.nunique())},
    {"statistic": "n_target_comparisons",
     "value": int(len(pairs)) if len(pairs) else 0},
    {"statistic": "n_comparisons_incomplete",
     "value": int((pairs.status == "incomplete").sum()) if len(pairs) else 0},
    {"statistic": "n_comparisons_not_independent", "value": n_dep},
    {"statistic": "n_comparisons_counted", "value": n_cmp},
    {"statistic": "n_comparisons_with_two_cq_values", "value": n_informative},
    {"statistic": "n_calls_agreeing", "value": n_agree},
    # Agreement over pairs that both returned no amplification is 1 by
    # construction. It is only a reproducibility statistic where at least one
    # comparison has two independent Cq values to disagree about.
    {"statistic": "call_agreement",
     "value": (round(n_agree / n_cmp, 4) if n_cmp and n_informative
               else "not estimable")},
    {"statistic": "mean_abs_delta_cq",
     "value": (round(mad, 3) if pd.notna(mad) else "not estimable")},
    {"statistic": "max_abs_delta_cq",
     "value": (round(float(done.delta_cq.abs().max()), 3)
               if n_informative else "not estimable")},
]
if n_dep:
    summary.append({
        "statistic": "replicate_arm_verdict",
        "value": (f"{n_dep} of {len(pairs)} comparisons carry identical Cq "
                  f"and identical Cq-confidence across the two plates, so "
                  f"the second record is a copy rather than a second "
                  f"measurement. The inter-plate replicate arm does not "
                  f"estimate plate-to-plate reproducibility in this dataset "
                  f"and must not be cited as evidence of it.")})
elif n_informative:
    summary.append({
        "statistic": "replicate_arm_verdict",
        "value": (f"{n_informative} comparison"
                  f"{'' if n_informative == 1 else 's'} "
                  f"ha{'s' if n_informative == 1 else 've'} two independent "
                  f"Cq values; the agreement and Cq difference above are "
                  f"interpretable as plate-to-plate reproducibility.")})
else:
    summary.append({
        "statistic": "replicate_arm_verdict",
        "value": ("every counted comparison returned no amplification on "
                  "both plates, so the pairs agree by construction and "
                  "carry no information about reproducibility.")})

# ---- first-round comparison, when that round's calls are supplied ----
if a.round1_calls and a.round1_calls.exists():
    r1 = pd.read_csv(a.round1_calls)
    key = "sample_code"
    if key not in r1.columns:
        print(f"ERROR --round1-calls has no '{key}' column; columns: "
              f"{list(r1.columns)}", file=sys.stderr)
        sys.exit(4)
    r1[key] = r1[key].astype(str).str.replace("_", "", regex=False).str.upper()
    now = (calls[calls.role.eq("sample")]
           .assign(sample_code=calls.sample_code.astype(str).str.upper())
           .pivot_table(index="sample_code", columns="target", values="call",
                        aggfunc="first"))
    n_tot = n_ag = 0
    for t in targets:
        col = f"{t}_call"
        if col not in r1.columns or t not in now.columns:
            continue
        j = now[[t]].join(r1.set_index(key)[[col]], how="inner").dropna()
        n_tot += len(j)
        n_ag += int((j[t] == j[col]).sum())
    summary += [
        {"statistic": "round1_n_comparisons", "value": n_tot},
        {"statistic": "round1_call_agreement",
         "value": round(n_ag / n_tot, 4) if n_tot else np.nan},
    ]
else:
    n_known = int(asg.first_run_plate.astype(str)
                  .ne("not run in first round").sum())
    summary += [
        {"statistic": "round1_n_specimens_with_first_run_plate", "value": n_known},
        {"statistic": "round1_call_agreement",
         "value": "not compared -- supply --round1-calls"},
    ]

pd.DataFrame(summary).to_csv(a.out_summary, index=False)
print(pd.DataFrame(summary).to_string(index=False))

if n_dep:
    print(f"\nWARNING {n_dep} of {len(pairs)} replicate comparisons are not "
          f"independent: the two plates carry identical Cq and identical "
          f"Cq-confidence, so the replicate record duplicates the primary "
          f"rather than measuring it again.", file=sys.stderr)
    print(pairs.loc[pairs.status == "not_independent",
                    ["sample_code", "target", "primary_plate", "primary_well",
                     "primary_cq", "replicate_plate", "replicate_well",
                     "replicate_cq"]].to_string(index=False), file=sys.stderr)
    print("  Check the plate exports for these wells against the instrument "
          "run files. Until they are resolved, the replicate arm does not "
          "support a reproducibility claim.", file=sys.stderr)

if n_cmp:
    print(f"\n{n_agree}/{n_cmp} counted target calls agree across the "
          f"inter-plate replicate pairs")
    bad = done[~done.call_agrees.astype(bool)]
    if len(bad):
        print("disagreeing comparisons:")
        print(bad[["sample_code", "target", "primary_plate", "primary_call",
                   "primary_cq", "replicate_plate", "replicate_call",
                   "replicate_cq", "delta_cq"]].to_string(index=False))
