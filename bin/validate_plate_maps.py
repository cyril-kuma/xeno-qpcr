#!/usr/bin/env python3
"""Confirm the plates that ran are the plates the design specified.

The randomisation fixes which specimen occupies which slot, and the slot fixes
three wells — one per target block. Nothing downstream can detect a plate set
up from the wrong row of the assignment: the calls would be internally
consistent and attached to the wrong specimen, and the zone, site and method a
specimen contributes to would all be wrong. This stage is the only place that
failure is visible.

Every specimen well in the merged table is matched against the assignment on
(plate, well, target). Four classes of disagreement abort the run:

  misplaced   the well holds a different specimen than the assignment expects
  unexpected  the export has a specimen well the assignment does not list
  absent      the assignment lists a well the export does not contain
  wrong_slot  a specimen appears on its assigned plate but not in its slot

Validation plates are exempt: they carry known-composition material rather
than field specimens, so the assignment does not describe them.

The validated well table is this stage's output, and the stages after it read
that file rather than the merged one. The dependency is what orders the
workflow, so the gate cannot be bypassed by running a later stage directly.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--wells", required=True, type=Path)
ap.add_argument("--plate-map", required=True, type=Path,
                help="expected plate x well x target x specimen, from "
                     "import_sample_sheet.py")
ap.add_argument("--out-wells", required=True, type=Path)
ap.add_argument("--out-report", required=True, type=Path)
ap.add_argument("--allow-plate-map-mismatch", action="store_true",
                help="record mismatches and continue. For diagnosing a "
                     "mis-exported plate; never for a reportable run")
a = ap.parse_args()

w = pd.read_csv(a.wells)
exp = pd.read_csv(a.plate_map)

routine = w[w.plate_class.eq("routine") & w.role.eq("sample")].copy()
routine["plate_id"] = routine.plate.astype(str).str.upper()
routine["well"] = routine.well.astype(str).str.upper()
exp["plate_id"] = exp.plate_id.astype(str).str.upper()
exp["well"] = exp.well.astype(str).str.upper()

# Only plates that were actually run can be checked; a design plate with no
# export is a missing run, which the report states rather than treats as a
# mismatch.
ran = set(routine.plate_id)
designed = set(exp.plate_id)
exp_ran = exp[exp.plate_id.isin(ran)]

m = exp_ran.merge(
    routine[["plate_id", "well", "target", "sample_code", "sample_name"]],
    on=["plate_id", "well", "target"], how="outer", indicator=True,
    suffixes=("_expected", "_observed"))

misplaced = m[(m._merge == "both")
              & (m.code.astype(str) != m.sample_code_observed.astype(str))]
unexpected = m[m._merge == "right_only"]
absent = m[m._merge == "left_only"]

# A specimen on the right plate but in the wrong slot: its code is somewhere
# on the plate, just not where the assignment put it. Worth naming separately,
# because the cause is a transcription slip rather than a swapped plate.
wrong_slot = misplaced[
    misplaced.apply(
        lambda r: r.code in set(
            routine.loc[routine.plate_id.eq(r.plate_id), "sample_code"]), axis=1)
] if len(misplaced) else misplaced

rows = [
    {"check": "plates_in_assignment", "n": len(designed)},
    {"check": "plates_ingested", "n": len(ran)},
    {"check": "plates_designed_not_run", "n": len(designed - ran),
     "detail": ",".join(sorted(designed - ran))},
    {"check": "plates_run_not_in_assignment", "n": len(ran - designed),
     "detail": ",".join(sorted(ran - designed))},
    {"check": "specimen_wells_expected", "n": len(exp_ran)},
    {"check": "specimen_wells_observed", "n": len(routine)},
    {"check": "wells_matched", "n": int((m._merge == "both").sum())},
    {"check": "wells_misplaced", "n": len(misplaced)},
    {"check": "wells_wrong_slot_same_plate", "n": len(wrong_slot)},
    {"check": "wells_unexpected", "n": len(unexpected)},
    {"check": "wells_absent", "n": len(absent)},
]
rep = pd.DataFrame(rows)
rep["detail"] = rep.get("detail", pd.Series(dtype=str)).fillna("")
rep.to_csv(a.out_report, index=False)
print(rep.to_string(index=False))

problems = {"misplaced": misplaced, "unexpected": unexpected, "absent": absent}
n_bad = sum(len(v) for v in problems.values())

if n_bad:
    print(f"\n{'ERROR' if not a.allow_plate_map_mismatch else 'WARNING'} "
          f"the plates that ran do not match the assignment "
          f"({n_bad} wells).", file=sys.stderr)
    for name, d in problems.items():
        if not len(d):
            continue
        cols = ["plate_id", "well", "target", "code", "sample_code_observed"]
        print(f"\n  {name} ({len(d)}):", file=sys.stderr)
        print(d[cols].head(25).to_string(index=False), file=sys.stderr)
        if len(d) > 25:
            print(f"  ... and {len(d) - 25} more", file=sys.stderr)
    if len(wrong_slot):
        print(f"\n  {len(wrong_slot)} of the misplaced wells hold a specimen "
              f"that is on its assigned plate but in the wrong slot, which "
              f"points at a transcription slip rather than a swapped plate.",
              file=sys.stderr)
    if not a.allow_plate_map_mismatch:
        print("\n  Resolve against the sample sheet's assignment tab and "
              "the laboratory run log before analysing. Calls from a mis-set "
              "plate are internally consistent and attached to the wrong "
              "specimen, so no later stage can detect this.", file=sys.stderr)
        sys.exit(6)

w.to_csv(a.out_wells, index=False)
n_val = int(w.plate_class.eq("validation").sum())
print(f"\nvalidated {len(routine)} specimen wells on {len(ran)} routine plates"
      + (f"; {n_val} wells on validation plates carried through unchecked"
         if n_val else ""))
