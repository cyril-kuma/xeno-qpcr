#!/usr/bin/env python3
"""Ingest one QuantStudio export into a tidy well table.

One plate per invocation so Nextflow can fan out across plates.

INGEST INTEGRITY
----------------
An earlier analysis of this dataset silently lost the NED/Po channel on 42 of
60 plates: wells the instrument had called `Amp` arrived as `No Amp` with a
null Cq and a Cq-confidence of exactly 0.000, which pushed reported P. ovale
prevalence from ~60% down to 4.3%. Nothing in that pipeline noticed.

This script therefore re-reads the raw sheet and asserts, per target, that the
ingested amplification-status counts equal the raw counts. Any drift aborts
the run rather than propagating. `--allow-ingest-drift` downgrades the abort
to a warning for debugging only; never use it for a reportable run.

ROLE INTEGRITY
--------------
Roles come from `Sample Name` rather than `Task`, because `Task == "NTC"`
covers both no-template controls and extraction blanks (see well_roles.py).
Every well must resolve to a role, and that role must be one the well's `Task`
can carry; a disagreement aborts the run.

PLATE CLASS
-----------
Routine analytical plates and validation plates are ingested by the same code
and separated by `plate_class`, taken from the plate identifier. Validation
plates carry known-composition specimens and characterise the assay; they are
excluded from prevalence and feed the cross-reactivity assessment instead.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

from well_roles import SAMPLE, assign_role, control_target

RENAME = {
    "Well": "well_index", "Well Position": "well", "Omit": "omit",
    "Sample Name": "sample_name", "Target Name": "target", "Task": "task",
    "Reporter": "reporter", "CT": "cq", "Quantity": "quantity",
    "Cq Conf": "cq_conf", "Amp Status": "amp_status",
    "Y-Intercept": "y_intercept", "R(superscript 2)": "r2",
    "Slope": "slope", "Efficiency": "efficiency",
}
NUMERIC = ["cq", "quantity", "cq_conf", "y_intercept", "r2", "slope", "efficiency"]

#: Plate identifier in a file name: a letter prefix, a run number, and an
#: optional repeat suffix, so Plas_18S_qPCR_R07-2.xlsx is plate R07-2 rather
#: than plate 2.
PLATE_ID_RE = re.compile(r"([A-Za-z]+\d+(?:-\d+)?)$")


def plate_id_from_path(path: Path) -> str:
    """Plate identifier from a file name, keeping any repeat suffix."""
    m = PLATE_ID_RE.search(path.stem)
    if not m:
        raise SystemExit(
            f"ERROR {path.name}: no plate identifier at the end of the file "
            f"name. Expected a letter prefix, a number and an optional "
            f"repeat suffix, e.g. Plas_18S_qPCR_R07.xlsx or "
            f"Plas_18S_qPCR_R07-2.xlsx."
        )
    return m.group(1).upper()


def read_sheet(path: Path, sheet: str, header):
    """Read with calamine first: 57/60 QuantStudio exports carry a malformed
    workbookView attribute that openpyxl rejects outright."""
    last = None
    for engine in ("calamine", "openpyxl"):
        try:
            return pd.read_excel(path, sheet_name=sheet, header=header, engine=engine)
        except Exception as exc:  # noqa: BLE001
            last = exc
    raise RuntimeError(f"{path.name}: unreadable by calamine or openpyxl ({last})")


def find_header(path: Path, sheet: str, fallback: int) -> int:
    raw = read_sheet(path, sheet, None)
    for i in range(min(len(raw), 150)):
        row = raw.iloc[i].astype(str).tolist()
        if "Well" in row and "Target Name" in row:
            return i
    return fallback


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--plate", required=True, type=Path)
    ap.add_argument("--sheet", default="Results")
    ap.add_argument("--header-row", type=int, default=46)
    ap.add_argument("--validation-plate-pattern", default=r"^V",
                    help="regex on the plate identifier marking a validation "
                         "plate; matching plates get plate_class=validation "
                         "and are excluded from prevalence")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--audit", required=True, type=Path)
    ap.add_argument("--allow-ingest-drift", action="store_true")
    args = ap.parse_args()

    plate_id = plate_id_from_path(args.plate)
    plate_class = ("validation"
                   if re.search(args.validation_plate_pattern, plate_id, re.I)
                   else "routine")
    hdr = find_header(args.plate, args.sheet, args.header_row)
    raw = read_sheet(args.plate, args.sheet, hdr)
    raw = raw[raw["Well"].notna()].copy()

    # ---- reference counts straight off the sheet, before any processing ----
    ref = (raw.assign(_t=raw["Target Name"].astype(str).str.strip(),
                      _a=raw["Amp Status"].astype(str).str.strip())
              .groupby(["_t", "_a"]).size().to_dict())

    df = raw.rename(columns=RENAME)
    df = df[[c for c in RENAME.values() if c in df.columns]].copy()
    for c in NUMERIC:
        if c in df.columns:
            df[c] = pd.to_numeric(
                df[c].replace({"Undetermined": None, "undetermined": None}),
                errors="coerce")
    for c in ("target", "task", "amp_status", "sample_name", "reporter"):
        if c in df.columns:
            df[c] = df[c].astype(str).str.strip()
    df["plate"] = plate_id
    df["plate_class"] = plate_class

    # ---- roles from Sample Name, cross-checked against Task ----
    assigned = [assign_role(n, t) for n, t in zip(df["sample_name"], df["task"])]
    df["role"] = [r for r, _ in assigned]
    consistency = [c for _, c in assigned]
    df["control_target"] = [control_target(n, r) for n, r
                            in zip(df["sample_name"], df["role"])]

    role_problems = sorted({c for c in consistency if c != "ok"})
    if role_problems:
        bad = df.loc[[c != "ok" for c in consistency],
                     ["well", "sample_name", "task", "target"]]
        print(f"ERROR [{plate_id}] ROLE INTEGRITY: Sample Name and Task "
              f"disagree, or a well resolves to no role.", file=sys.stderr)
        for p in role_problems:
            print(f"  - {p}", file=sys.stderr)
        print(bad.to_string(index=False), file=sys.stderr)
        print("  Roles are assigned from Sample Name (see bin/well_roles.py). "
              "Correct the export's Sample Name or extend the role patterns; "
              "do not relax the check.", file=sys.stderr)
        return 5

    # canonical join key: metadata uses NS_141, the plate export uses NS141
    df["sample_code"] = (df["sample_name"].str.replace("_", "", regex=False)
                         .str.upper())
    df.loc[df["role"].ne(SAMPLE), "sample_code"] = pd.NA

    # ---- integrity assertion ----
    got = df.groupby(["target", "amp_status"]).size().to_dict()
    drift = {f"{t}|{a}": {"raw": int(n), "ingested": int(got.get((t, a), 0))}
             for (t, a), n in ref.items() if int(got.get((t, a), 0)) != int(n)}

    audit = {
        "plate": plate_id, "plate_class": plate_class,
        "file": args.plate.name, "header_row": int(hdr),
        "n_wells": int(len(df)),
        "targets": sorted(df["target"].unique().tolist()),
        "n_by_role": {k: int(v) for k, v in df["role"].value_counts().items()},
        "amp_status_by_target": {f"{t}|{a}": int(n) for (t, a), n in ref.items()},
        "ingest_drift": drift,
        "status": "ok" if not drift else "DRIFT",
    }
    args.audit.write_text(json.dumps(audit, indent=2))
    df.to_csv(args.out, index=False)

    if drift:
        msg = (f"[{plate_id}] INGEST DRIFT: raw and ingested amplification-status "
               f"counts disagree: {drift}")
        if not args.allow_ingest_drift:
            print("ERROR " + msg, file=sys.stderr)
            return 2
        print("WARNING " + msg, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
