#!/usr/bin/env python3
"""Prevalence tables with Wilson intervals plus a run report."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportion_confint

ap = argparse.ArgumentParser()
ap.add_argument("--mosquito", required=True, type=Path)
ap.add_argument("--sensitivity", required=True, type=Path)
ap.add_argument("--target-qc", required=True, type=Path)
ap.add_argument("--control-qc", required=True, type=Path)
ap.add_argument("--diagnostic", required=True, type=Path)
ap.add_argument("--allocation-balance", required=True, type=Path)
ap.add_argument("--ingest-audit", required=True, type=Path)
ap.add_argument("--plate-map-report", required=True, type=Path)
ap.add_argument("--replicate-summary", required=True, type=Path)
ap.add_argument("--ipc-summary", required=True, type=Path)
ap.add_argument("--cross-reactivity", required=True, type=Path)
ap.add_argument("--reconciliation", required=True, type=Path)
ap.add_argument("--version", default="")
ap.add_argument("--zone-col", default="bioclimatic_zone")
ap.add_argument("--out-report", required=True, type=Path)
ap.add_argument("--out-prevalence", required=True, type=Path)
a = ap.parse_args()

m = pd.read_csv(a.mosquito)
sens = pd.read_csv(a.sensitivity)
qc = pd.read_csv(a.target_qc)
ctl = pd.read_csv(a.control_qc)
diag = pd.read_csv(a.diagnostic)
bal = pd.read_csv(a.allocation_balance)
aud = pd.read_csv(a.ingest_audit)
pmap = pd.read_csv(a.plate_map_report)
rep = pd.read_csv(a.replicate_summary)
ipc = pd.read_csv(a.ipc_summary)
crx = pd.read_csv(a.cross_reactivity)
rec = pd.read_csv(a.reconciliation)
targets = sens.target.unique().tolist()


def ci(k, n):
    if n == 0:
        return (np.nan, np.nan)
    return proportion_confint(k, n, alpha=0.05, method="wilson")


rows = []


def add(zone, outcome, elig, det):
    n, k = int(elig.sum()), int((elig & det).sum())
    lo, hi = ci(k, n)
    rows.append({"zone": zone, "outcome": outcome, "n_eligible": n,
                 "n_detected": k,
                 "prevalence": round(k / n, 4) if n else np.nan,
                 "wilson_lo": round(lo, 4), "wilson_hi": round(hi, 4)})


groups = [("Overall", m)]
if a.zone_col in m.columns:
    groups += [(z, m[m[a.zone_col] == z])
               for z in sorted(m[a.zone_col].dropna().unique())]
for zone, d in groups:
    add(zone, "any_plasmodium", d.eligible_any, d.any_plasmodium_detected)
    for t in targets:
        add(zone, t, d[f"{t}_eligible"], d[f"{t}_detected"])
    add(zone, "mixed", d.eligible_any, d.mixed_detected)
prev = pd.DataFrame(rows)
prev.to_csv(a.out_prevalence, index=False)


def stat(df, name, default="n/a"):
    hit = df.loc[df.statistic == name, "value"]
    return hit.iloc[0] if len(hit) else default


def get(df, col, name, default="n/a"):
    hit = df.loc[df[col] == name, "value" if "value" in df.columns else "n"]
    return hit.iloc[0] if len(hit) else default


n_routine = int((aud.plate_class == "routine").sum())
n_val = int((aud.plate_class == "validation").sum())
flag = str(stat(diag, "batch_effect_flag"))
note = str(stat(diag, "interpretation"))
crx_verdict = str(get(crx, "check", "validation_plates", "0"))

L = [f"# xeno-qpcr run report{(' ' + a.version) if a.version else ''}", "",
     "## Run integrity", "",
     f"- plates ingested: {len(aud)} "
     f"({n_routine} routine, {n_val} validation); "
     f"ingest drift on {int((aud.status != 'ok').sum())}",
     f"- specimen wells matched to the assignment: "
     f"{get(pmap, 'check', 'wells_matched')} of "
     f"{get(pmap, 'check', 'specimen_wells_expected')} expected "
     f"(misplaced {get(pmap, 'check', 'wells_misplaced')}, "
     f"unexpected {get(pmap, 'check', 'wells_unexpected')}, "
     f"absent {get(pmap, 'check', 'wells_absent')})",
     f"- plate x target standard curves passing MIQE QC: "
     f"{int(qc.curve_ok.sum())}/{len(qc)}",
     f"- specimens analysed: {len(m)}", "",
     "### Specimens", "",
     "The sample sheet is reconciled in both directions rather than "
     "intersected, so the denominator chain is explicit: specimens run, "
     "specimens with a record, specimens analysed.", "",
     rec[["stage", "n"]].to_markdown(index=False), ""]
n_un = int(rec.loc[rec.stage.eq("assigned but not described"), "n"].iloc[0]) \
    if (rec.stage == "assigned but not described").any() else 0
if n_un:
    ids = rec.loc[rec.stage.eq("assigned but not described"), "detail"].iloc[0]
    L += [f"> {n_un} specimens were run but have no record in the specimens "
          f"tab, so they are excluded from every denominator: "
          f"`{ids}`.", ""]
L += []

# ---- controls, one line per class ----
L += ["### Controls", "",
      "Each control class is gated separately; the failure each implicates "
      "differs, and so does the corrective action.", ""]
cl = (ctl.assign(fail=ctl.verdict.str.startswith("FAIL"),
                 warn=ctl.verdict.str.startswith("WARN"))
         .groupby("role")
         .agg(n_wells=("verdict", "size"), n_fail=("fail", "sum"),
              n_warn=("warn", "sum"))
         .reset_index())
L += [cl.to_markdown(index=False), ""]

L += ["### Cross-reactivity", "",
      crx[["check", "n", "detail"]].to_markdown(index=False), ""]

L += ["## Prevalence", "", prev.to_markdown(index=False), "",
      "## Cross-amplification sensitivity", "",
      "The guardband is a sensitivity arm and never alters the primary call.",
      "", sens.to_markdown(index=False), ""]

L += ["## Allocation balance", "",
      "Plate against each stratifying variable. The randomisation exists so "
      "that plate is not confounded with site, zone, method or year; this is "
      "the check that it worked.", "",
      bal.to_markdown(index=False), ""]

L += ["## Residual plate structure", "", diag.to_markdown(index=False), ""]

if flag == "CONFOUNDED":
    L += ["> **The allocation is imbalanced.** Plate cannot be separated from "
          "the flagged variable, so no plate-level result is interpretable "
          "until the allocation is resolved.", ""]
elif flag == "PLATE_EFFECT":
    L += ["> **Detection is plate-clustered, and the allocation is balanced.**",
          f"> {note}.", ""]
else:
    L += [f"> {note[0].upper() + note[1:]}.", ""]

L += ["## Inter-plate replicate concordance", "",
      "Six specimens were aliquoted onto a second plate. These pairs are the "
      "only measurement in the design that separates plate-to-plate technical "
      "variation from specimen-to-specimen biological variation.", "",
      rep.to_markdown(index=False), ""]

L += ["## Internal positive control", "",
      "Without an IPC a well with no Cq cannot be distinguished from a failed "
      "extraction or an inhibited reaction (assay design §10). The count of "
      "specimens returning no signal at every target is reported here so the "
      "effect can be bounded.", "",
      ipc.to_markdown(index=False), ""]

L += ["---", "",
      "Because extracts are whole-mosquito, a positive means *Plasmodium* DNA "
      "was present in the mosquito. It does not distinguish parasites ingested "
      "with the blood meal from an established infection, and it is not a "
      "sporozoite rate.", ""]

Path(a.out_report).write_text("\n".join(L))
print(f"wrote {a.out_report}")
print(f"  allocation: {'balanced' if stat(diag, 'allocation_decorrelated') else 'IMBALANCED'}"
      f" | plate structure: {flag} | cross-reactivity plates: {crx_verdict}")
