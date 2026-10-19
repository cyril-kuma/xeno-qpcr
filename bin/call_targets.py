#!/usr/bin/env python3
"""Well-level calls and mosquito-level collapse.

PRIMARY RULE
------------
  amp_status == "Inconclusive"          -> invalid
  amp_status == "No Amp" or Cq is null  -> negative
  Cq <= cq_positive_max (38)            -> positive
        ...except targets in cqconf_targets (Pm), which additionally require
           Cq Conf >= cqconf_min (0.50); otherwise trace
  cq_positive_max < Cq <= cq_trace_max  -> trace
  Cq > cq_trace_max (40)                -> negative

Trace is NOT positive. It is retained as its own state so the evaluable
denominator stays honest.

ELIGIBILITY
-----------
  per-target eligible  = call != invalid
  eligible_any         = (>=1 target positive) OR (all targets valid)
  eligible_all         = all targets valid            (mixed-detection denominator)

CROSS-AMPLIFICATION GUARDBAND
-----------------------------
Sensitivity arm only; it never touches the primary call. A non-dominant
target positive on a mosquito whose dominant target (Pf) is also positive,
and whose Cq is >= dominant Cq + guardband, is demoted. Applied to every
non-dominant target.

INTERNAL POSITIVE CONTROL
-------------------------
Without an IPC a well with no Cq cannot be distinguished from a failed
extraction or an inhibited reaction, so it cannot be called a true negative
(assay design §10). The IPC was not incorporated during the study period, so
`--ipc-enabled` is off and the output instead carries
`all_targets_undetermined`, the count assay design §10.6 requires to be
reported so readers can bound the effect. When an IPC is added, enabling it
routes a specimen whose IPC fails to `indeterminate` at every target rather
than to negative.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

#: Wells whose amplification status means the instrument saw nothing, as
#: opposed to seeing something above the trace band.
NO_SIGNAL = {"No Amp"}


def normalise_key(series: pd.Series, pattern: str | None = None) -> pd.Series:
    """Canonical join key: optional regex extraction, then strip and upper-case.

    Plate exports carry the bare Sample Name (NS_003). Metadata files sometimes
    carry a composite (NS_003_NIH00124L). The pattern pulls the joinable part
    out before normalisation.
    """
    v = series.astype(str).str.strip()
    if pattern:
        v = v.str.extract(pattern, expand=False).fillna(series.astype(str))
    return v.str.replace("_", "", regex=False).str.upper()


DEFAULT_PATTERNS = [
    (r"^([A-Za-z]+_?\d+)", "leading letters + digits, e.g. NS_003 from "
                            "NS_003_NIH00124L"),
]


def suggest_key(meta: pd.DataFrame, plate_codes: set, tried: str) -> str:
    """Scan the metadata for a column, or a pattern, that would actually join."""
    best = []
    for col in meta.columns:
        if col == "sample_code":
            continue
        raw = normalise_key(meta[col])
        n = len(plate_codes & set(raw))
        if n:
            best.append((n, col, None))
        for pat, _desc in DEFAULT_PATTERNS:
            ext = normalise_key(meta[col], pat)
            n2 = len(plate_codes & set(ext))
            if n2 > n:
                best.append((n2, col, pat))
    if not best:
        return ("  No column in the metadata joins to the plate codes, with or "
                "without pattern extraction. The two files describe different "
                "specimens.")
    best.sort(reverse=True)
    n, col, pat = best[0]
    if pat is None:
        return (f"  '{tried}' does not join, but '{col}' matches {n} plate "
                f"codes. Re-run with --metadata-key {col}\n"
                f"  (Nextflow: --metadata_key {col})")
    return (f"  '{tried}' does not join as-is, but extracting with {pat!r} "
            f"matches {n} plate codes -- the column looks like a composite "
            f"identifier.\n"
            f"  Re-run with: --metadata-key {col} "
            f"--metadata-key-pattern '{pat}'\n"
            f"  (Nextflow: --metadata_key {col} "
            f"--metadata_key_pattern '{pat}')")


def _sample_codes(df, n: int = 5) -> list:
    """A few example identifiers, NaN-safe and sortable."""
    vals = df["sample_code"].dropna().astype(str).unique().tolist()
    return sorted(vals)[:n]


def die(msg: str) -> None:
    """Abort with an actionable message rather than a stack trace."""
    print(f"ERROR {msg}", file=sys.stderr)
    sys.exit(4)


def base_call(amp_status, cq, target, cq_conf, p) -> str:
    if amp_status == "Inconclusive":
        return "invalid"
    if amp_status in NO_SIGNAL or pd.isna(cq):
        return "negative"
    if amp_status != "Amp":
        return "invalid"
    if cq <= p["cq_positive_max"]:
        if target in p["cqconf_targets"]:
            if pd.isna(cq_conf) or cq_conf < p["cqconf_min"]:
                return "trace"
        return "positive"
    if cq <= p["cq_trace_max"]:
        return "trace"
    return "negative"


def apply_codetection_exclusion(wide, targets, dominant, exclude):
    """Primary specificity rule: a target in `exclude` cannot be called
    positive in a mosquito where the dominant target is positive.

    Rationale. The assay follows the Shokoples design: species-specific
    forward primers with a shared conserved reverse primer for Pm and Po,
    leaving discrimination to the forward primer and a short MGB probe.
    Empirically 92.6% of Po positives co-occur with Pf at a median delay of
    +3.1 cycles, the signature of the Po reaction amplifying Pf template. Any
    co-detection is therefore treated as cross-amplification rather than as a
    mixed infection, and reported as such.

    This is stricter than a guardband of 0, which retains co-detections whose
    Cq happens to fall below the dominant target's.
    """
    out = wide.copy()
    dom_pos = out[f"{dominant}_call"].eq("positive")
    for t in targets:
        if t == dominant or t not in exclude:
            continue
        hit = out[f"{t}_call"].eq("positive") & dom_pos
        out.loc[hit, f"{t}_call"] = "excluded_codetection"
    return out


def apply_guardband(wide: pd.DataFrame, targets, dominant: str,
                    guardband: float) -> pd.DataFrame:
    """Demote every non-dominant positive delayed by >= guardband cycles."""
    out = wide.copy()
    dom_pos = out[f"{dominant}_call"].eq("positive")
    dom_cq = out[f"{dominant}_cq"]
    for t in targets:
        if t == dominant:
            out[f"{t}_call_gb"] = out[f"{t}_call"]
            continue
        delta = out[f"{t}_cq"] - dom_cq
        out[f"{t}_delta_cq_vs_{dominant.lower()}"] = delta
        demote = (out[f"{t}_call"].eq("positive") & dom_pos
                  & delta.notna() & (delta >= guardband))
        out[f"{t}_call_gb"] = np.where(demote, "demoted_crossamp", out[f"{t}_call"])
    return out


def apply_ipc(wide: pd.DataFrame, targets, ipc_cq_max: float) -> pd.DataFrame:
    """Route specimens whose IPC failed away from negative.

    A specimen with no IPC amplification, or IPC amplification later than
    `ipc_cq_max`, has not demonstrated that its extraction worked and that the
    reaction was uninhibited. Its target calls become `indeterminate` and it
    drops out of every denominator; it is never reported as negative.
    Positive calls stand: a target that amplified did so regardless of what
    the IPC shows.
    """
    out = wide.copy()
    ok = out["ipc_cq"].notna() & (out["ipc_cq"] <= ipc_cq_max)
    out["ipc_pass"] = ok
    failed = ~ok
    for t in targets:
        hit = failed & out[f"{t}_call"].ne("positive")
        out.loc[hit, f"{t}_call"] = "indeterminate"
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--wells", required=True, type=Path)
    ap.add_argument("--metadata", type=Path)
    ap.add_argument("--targets", default="Pf,Pm,Po")
    ap.add_argument("--dominant", default="Pf")
    ap.add_argument("--cq-positive-max", type=float, default=38.0)
    ap.add_argument("--cq-trace-max", type=float, default=40.0)
    ap.add_argument("--cqconf-targets", default="Pm",
                    help="comma-separated targets that must also clear the "
                         "Cq-confidence floor; empty string disables")
    ap.add_argument("--cqconf-min", type=float, default=0.50)
    ap.add_argument("--guardband-grid", default="0,1,2,3,5")
    ap.add_argument("--codetection-exclude", default="Po",
                    help="targets that cannot be positive when the dominant "
                         "target is positive; empty string disables")
    ap.add_argument("--cross-reactivity", type=Path,
                    help="verdict from the validation plates, recorded "
                         "alongside the calls as the evidence state the "
                         "co-detection rule rests on")
    ap.add_argument("--ipc-enabled", action="store_true",
                    help="route specimens failing the internal positive "
                         "control to indeterminate rather than negative")
    ap.add_argument("--ipc-target", default="IPC",
                    help="Target Name of the internal positive control")
    ap.add_argument("--ipc-cq-max", type=float, default=32.0,
                    help="IPC Cq above which a specimen has not demonstrated "
                         "a successful, uninhibited reaction")
    ap.add_argument("--restrict-to-metadata", action="store_true")
    ap.add_argument("--metadata-key", default="sample_id",
                    help="metadata column holding the plate Sample Name")
    ap.add_argument("--metadata-key-pattern", default=None,
                    help="regex applied to --metadata-key; group 1 becomes the "
                         r"join key. Use when the metadata column is a composite, "
                         r"e.g. '^([A-Za-z]+_?\d+)' turns NS_003_NIH00124L into "
                         "NS_003")
    ap.add_argument("--out-calls", required=True, type=Path)
    ap.add_argument("--out-mosquito", required=True, type=Path)
    ap.add_argument("--out-sensitivity", required=True, type=Path)
    ap.add_argument("--out-ipc-summary", type=Path,
                    help="counts assay design §10.6 requires while the IPC is "
                         "unavailable")
    args = ap.parse_args()

    targets = [t.strip() for t in args.targets.split(",") if t.strip()]
    p = {"cq_positive_max": args.cq_positive_max,
         "cq_trace_max": args.cq_trace_max,
         "cqconf_min": args.cqconf_min,
         "cqconf_targets": {t.strip() for t in args.cqconf_targets.split(",")
                            if t.strip()}}

    w = pd.read_csv(args.wells)
    print(f"[1] wells: {len(w)} rows, columns={list(w.columns)}")
    for col in ("role", "target", "sample_code", "amp_status", "cq"):
        if col not in w.columns:
            die(f"the well table has no '{col}' column. Columns present: "
                f"{list(w.columns)}. This file is produced by "
                f"VALIDATE_PLATE_MAPS; check its work directory.")
    print(f"[2] role values: {w['role'].value_counts().to_dict()}")
    print(f"[3] target values: {w['target'].value_counts().to_dict()}")

    if "plate_class" in w.columns:
        n_val = int(w["plate_class"].eq("validation").sum())
        if n_val:
            print(f"[3a] excluding {n_val} wells on validation plates; "
                  f"prevalence is estimated from routine plates only")
        w = w[w["plate_class"].eq("routine")]

    s = w[w["role"] == "sample"].copy()
    if s.empty:
        die("no wells with role=='sample'. Roles are assigned from Sample "
            f"Name (see bin/well_roles.py). role values found: "
            f"{w['role'].value_counts().to_dict()}")
    print(f"[4] sample wells: {len(s)} across {s['sample_code'].nunique()} specimens")

    seen = set(s["target"].dropna().unique())
    unknown = seen - set(targets) - {args.ipc_target}
    absent = [t for t in targets if t not in seen]
    if absent:
        die(f"the well table contains no rows for target(s) {absent}. "
            f"Targets present: {sorted(seen)}. Either --targets does not match "
            f"the Target Name values in the plate exports, or those wells were "
            f"lost at ingest.")
    if unknown:
        print(f"[!] target names present but not requested: {sorted(unknown)}")

    if args.metadata and args.metadata.exists():
        meta = pd.read_csv(args.metadata)
        key = (args.metadata_key if args.metadata_key in meta.columns
               else meta.columns[0])
        meta["sample_code"] = normalise_key(meta[key], args.metadata_key_pattern)
        plate_codes = set(s["sample_code"].dropna())
        s["metadata_linked"] = s["sample_code"].isin(set(meta["sample_code"]))
        n_link = int(s["metadata_linked"].sum())
        print(f"[5] metadata: {len(meta)} rows, joining on '{key}'"
              + (f" via pattern {args.metadata_key_pattern!r}"
                 if args.metadata_key_pattern else "")
              + f"; {n_link}/{len(s)} sample wells linked")
        if n_link == 0:
            die("no sample well matched the metadata.\n"
                f"  Example plate codes:    {_sample_codes(s)}\n"
                f"  Example metadata codes: {_sample_codes(meta)}\n"
                + suggest_key(meta, plate_codes, key))
        if args.restrict_to_metadata:
            s = s[s["metadata_linked"]]
            print(f"[6] restricted to metadata: {len(s)} wells, "
                  f"{s['sample_code'].nunique()} specimens")
    else:
        meta = None
        s["metadata_linked"] = True
        print("[5] no metadata file supplied; all wells retained")

    if s.empty:
        die("no sample wells remain after filtering. Re-run without "
            "--restrict-to-metadata to confirm the metadata join is the cause.")

    ipc_wells = s[s["target"].eq(args.ipc_target)].copy()
    s = s[s["target"].ne(args.ipc_target)]

    s["call"] = [base_call(a, c, t, q, p) for a, c, t, q in
                 zip(s["amp_status"], s["cq"], s["target"],
                     s.get("cq_conf", pd.Series(np.nan, index=s.index)))]
    s["no_signal"] = s["amp_status"].isin(NO_SIGNAL) | s["cq"].isna()
    s.to_csv(args.out_calls, index=False)

    # The wide table is built explicitly, one target at a time. Unpacking a
    # pivot_table MultiIndex gives column names whose order is not guaranteed
    # across pandas versions, and a target with no rows is dropped silently.
    wide = pd.DataFrame({"sample_code": s["sample_code"].drop_duplicates()
                         .sort_values().to_numpy()})
    wide = wide.set_index("sample_code")
    for t in targets:
        sub = (s[s["target"].eq(t)]
               .drop_duplicates("sample_code")
               .set_index("sample_code"))
        wide[f"{t}_call"] = sub["call"]
        wide[f"{t}_cq"] = sub["cq"]
        wide[f"{t}_no_signal"] = sub["no_signal"].reindex(wide.index).fillna(False)
    wide = wide.reset_index()

    missing = [t for t in targets if s["target"].eq(t).sum() == 0]
    if missing:
        print(f"WARNING no wells found for target(s) {missing}; their columns "
              f"are present but empty", file=sys.stderr)

    # ---- internal positive control ----
    ipc_note = "not incorporated during the study period (assay design §10)"
    if len(ipc_wells):
        ipc_cq = (ipc_wells.drop_duplicates("sample_code")
                  .set_index("sample_code")["cq"])
        wide["ipc_cq"] = wide["sample_code"].map(ipc_cq)
    else:
        wide["ipc_cq"] = np.nan

    if args.ipc_enabled:
        if not len(ipc_wells):
            die(f"--ipc-enabled was given but no wells carry Target Name "
                f"{args.ipc_target!r}. Targets present: {sorted(seen)}. Either "
                f"the IPC was not run on these plates, or --ipc-target does "
                f"not match the name the export uses.")
        wide = apply_ipc(wide, targets, args.ipc_cq_max)
        n_fail = int((~wide["ipc_pass"]).sum())
        ipc_note = (f"enabled; {n_fail} specimens failed the IPC "
                    f"(Cq > {args.ipc_cq_max} or absent) and are "
                    f"indeterminate rather than negative")
        print(f"[7] IPC: {int(wide['ipc_pass'].sum())}/{len(wide)} specimens "
              f"pass; {n_fail} routed to indeterminate")
    else:
        wide["ipc_pass"] = pd.NA

    # Specimens with no signal at any target. While the IPC is unavailable
    # these cannot be separated from failed extractions, and assay design
    # §10.6 requires the count to be reported.
    wide["all_targets_undetermined"] = wide[
        [f"{t}_no_signal" for t in targets]].all(axis=1)
    n_undet = int(wide["all_targets_undetermined"].sum())

    # ---- cross-reactivity evidence state ----
    verdict, verdict_detail = "not_assessed", (
        "no validation plate was supplied, so the absence of cross-reactivity "
        "with the uninfected Anopheles genome rests on the per-plate "
        "uninfected-host negative controls alone")
    if args.cross_reactivity and args.cross_reactivity.exists():
        cr = json.loads(args.cross_reactivity.read_text())
        verdict = cr.get("verdict", "not_assessed")
        verdict_detail = cr.get("detail", verdict_detail)
    wide["cross_reactivity_verdict"] = verdict
    print(f"[8] cross-reactivity verdict: {verdict} -- {verdict_detail}")

    # The guardband sensitivity arm asks what would happen if a Cq-delay
    # guardband were used INSTEAD of co-detection exclusion, so it is computed
    # on the pre-exclusion calls. Running it afterwards leaves it nothing to
    # act on and produces a flat column.
    pre_exclusion = wide.copy()

    exclude = {t.strip() for t in args.codetection_exclude.split(",") if t.strip()}
    if exclude:
        wide = apply_codetection_exclusion(wide, targets, args.dominant, exclude)

    for t in targets:
        wide[f"{t}_eligible"] = (wide[f"{t}_call"].notna()
                                 & ~wide[f"{t}_call"].isin(
                                     ["invalid", "indeterminate"]))
        wide[f"{t}_detected"] = wide[f"{t}_call"].eq("positive")

    n_pos = wide[[f"{t}_detected" for t in targets]].sum(axis=1)
    n_val_t = wide[[f"{t}_eligible" for t in targets]].sum(axis=1)
    wide["n_targets_positive"] = n_pos
    wide["n_targets_valid"] = n_val_t
    wide["any_plasmodium_detected"] = n_pos > 0
    wide["eligible_all"] = n_val_t == len(targets)
    wide["eligible_any"] = wide["any_plasmodium_detected"] | wide["eligible_all"]
    wide["mixed_detected"] = n_pos >= 2
    wide["detection_pattern"] = wide.apply(
        lambda r: "+".join(t for t in targets if r[f"{t}_detected"]) or "None",
        axis=1)

    if meta is not None:
        wide = wide.merge(meta, on="sample_code", how="left")
    wide.to_csv(args.out_mosquito, index=False)

    if args.out_ipc_summary:
        rows = [
            {"statistic": "ipc_status", "value": ipc_note},
            {"statistic": "n_specimens", "value": len(wide)},
            {"statistic": "n_all_targets_undetermined", "value": n_undet},
            {"statistic": "pct_all_targets_undetermined",
             "value": round(100 * n_undet / len(wide), 2) if len(wide) else np.nan},
        ]
        for t in targets:
            rows.append({"statistic": f"n_{t}_no_signal",
                         "value": int(wide[f"{t}_no_signal"].sum())})
        if args.ipc_enabled:
            rows.append({"statistic": "n_ipc_fail",
                         "value": int((~wide["ipc_pass"].astype(bool)).sum())})
        pd.DataFrame(rows).to_csv(args.out_ipc_summary, index=False)

    # ---- guardband sensitivity arm ----
    rows = []
    grid = [float(g) for g in args.guardband_grid.split(",") if g.strip()]
    for gb in grid:
        gbw = apply_guardband(pre_exclusion, targets, args.dominant, gb)
        for t in targets:
            elig = wide[f"{t}_eligible"]
            k = int((gbw[f"{t}_call_gb"].eq("positive") & elig).sum())
            n = int(elig.sum())
            rows.append({"guardband": gb, "target": t, "n_eligible": n,
                         "n_positive": k,
                         "prevalence": round(k / n, 4) if n else np.nan,
                         "n_demoted": int(gbw[f"{t}_call_gb"]
                                          .eq("demoted_crossamp").sum())})
    sens = pd.DataFrame(rows)
    sens.to_csv(args.out_sensitivity, index=False)

    # ---- no-op guard: a real guardband must move every non-dominant target --
    bad = []
    dominant_positive = int(wide[f"{args.dominant}_detected"].sum())
    for t in targets:
        if t == args.dominant or len(grid) < 2:
            continue
        # Nothing to demote is not a no-op: the guard is about the rule
        # failing to reach a target, not about an absence of co-detections.
        co = int((pre_exclusion[f"{t}_call"].eq("positive")
                  & pre_exclusion[f"{args.dominant}_call"].eq("positive")).sum())
        if co == 0 or dominant_positive == 0:
            print(f"[!] {t}: no co-detections with {args.dominant}; guardband "
                  f"sensitivity is uninformative for this target")
            continue
        if sens.loc[sens["target"].eq(t), "n_positive"].nunique() == 1:
            bad.append(t)
    if bad:
        print(f"ERROR guardband is a no-op for {bad}: identical counts at every "
              f"value in {grid}. The rule is not reaching these targets.",
              file=sys.stderr)
        return 3

    for t in targets:
        e = int(wide[f"{t}_eligible"].sum()); d = int(wide[f"{t}_detected"].sum())
        print(f"{t}: {d}/{e} ({100*d/e:.1f}%)" if e else f"{t}: no eligible")
    ea = int(wide["eligible_any"].sum())
    ap_ = int(wide["any_plasmodium_detected"].sum())
    print(f"any: {ap_}/{ea} ({100*ap_/ea:.1f}%)" if ea else "any: none")
    print(f"no signal at any target: {n_undet} specimens "
          f"(IPC {ipc_note.split(';')[0]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
