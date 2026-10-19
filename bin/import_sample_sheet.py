#!/usr/bin/env python3
"""Read the sample sheet: the specimens, and where each one was run.

`data/sample_sheet.xlsx` is the pipeline's only description of the experiment.
Two tabs:

  specimens    one row per specimen -- everything known about the mosquito,
               independent of where it was run
  assignment   one row per plate slot, from the randomisation. A slot fixes
               three wells, one per target block. Specimen attributes are not
               repeated here.

They are deliberately separate, because they are different entities at
different grain: 720 slots (714 specimens plus 6 inter-plate replicates)
against 701 specimens. Merging them would repeat specimen attributes on the
replicate rows and lose the distinction between specimens that were run and
specimens that are described.

The laboratory's own tracking workbook -- run log, QC checklist, reagent log,
manual results entry -- stays with the bench and is not read here. It is a
record of what was done, continues to be edited after an analysis is frozen,
and its manually transcribed results would be a second source of truth for
numbers the instrument exports already carry.

TWO CHECKS THIS STAGE MAKES
---------------------------
*Slot against wells.* `slot_id` and the three well columns are typed
independently but related by the Layout A block structure, so each row's well
triplet is checked against the mapping its slot implies. A mistyped well cell
is caught here rather than surfacing later as a plate that does not match its
own design.

*Specimens against assignment.* The two tabs are reconciled in both
directions and the counts reported, rather than silently intersected. A silent
inner join cannot distinguish a handful of known gaps from a mis-keyed column
that matches only half the sheet: both yield a clean run on whatever survived.
The reconciliation is what lets the analysis state its denominator chain --
N assigned, N described, N analysed -- and `--max-unmatched` stops a run where
the join has plainly broken.

Note that specimens with no record are **not** dropped here. They are
physically in wells on real plates, so the plate-map gate must see them;
removing them first would make their wells look unexpected. They leave the
analysis later, at CALL_TARGETS, under `restrict_to_metadata`.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

from randomise import slot_wells

TARGETS = ("Pf", "Pm", "Po")

SPECIMEN_REQUIRED = ["sample_code", "specimen_id", "bioclimatic_zone",
                     "collection_location", "collection_method",
                     "collection_date"]
ASSIGN_REQUIRED = ["plate_id", "run_order", "slot_id", "well_Pf", "well_Pm",
                   "well_Po", "sample_code", "stratum", "role"]

REPLICATE_RE = re.compile(r"INTER-PLATE REPLICATE \(of plate ([A-Za-z0-9-]+)\)",
                          re.I)


def read_tab(book: Path, tab: str, key: str) -> pd.DataFrame:
    """Read a tab whose header sits below a one-line note."""
    try:
        probe = pd.read_excel(book, tab, engine="calamine", header=None, nrows=4)
    except ValueError:
        raise SystemExit(
            f"ERROR {book.name} has no tab named {tab!r}. The sample sheet "
            f"needs a 'specimens' tab and an 'assignment' tab.")
    header = 0
    for i in range(min(4, len(probe))):
        if key in [str(v).strip() for v in probe.iloc[i].tolist()]:
            header = i
            break
    else:
        raise SystemExit(
            f"ERROR {book.name}!{tab}: no header row carrying {key!r} in the "
            f"first four rows.")
    return pd.read_excel(book, tab, engine="calamine", dtype=str, header=header)


def require(df: pd.DataFrame, cols: list[str], where: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR {where} is missing {missing}. Columns "
                         f"found: {list(df.columns)}")


def norm(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.replace("_", "", regex=False).str.upper()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-sheet", required=True, type=Path)
    ap.add_argument("--max-unmatched", type=int, default=20,
                    help="abort when more than this many assigned specimens "
                         "have no record in the specimens tab. A handful is a "
                         "documented gap; a large number means the join is "
                         "broken, and both look identical to a silent "
                         "intersection")
    ap.add_argument("--out-specimens", required=True, type=Path)
    ap.add_argument("--out-assignment", required=True, type=Path)
    ap.add_argument("--out-plate-map", required=True, type=Path)
    ap.add_argument("--out-reconciliation", required=True, type=Path)
    args = ap.parse_args()

    spec = read_tab(args.sample_sheet, "specimens", "sample_code")
    asg = read_tab(args.sample_sheet, "assignment", "plate_id")
    require(spec, SPECIMEN_REQUIRED, "the specimens tab")
    require(asg, ASSIGN_REQUIRED, "the assignment tab")

    spec = spec.dropna(subset=["sample_code"]).copy()
    asg = asg.dropna(subset=["plate_id", "slot_id"]).copy()

    spec["sample_code"] = spec.sample_code.str.strip()
    spec["code"] = norm(spec.sample_code)
    asg["plate_id"] = asg.plate_id.str.strip().str.upper()
    asg["sample_code"] = asg.sample_code.str.strip()
    asg["code"] = norm(asg.sample_code)
    for c in ("slot_id", "run_order"):
        asg[c] = pd.to_numeric(asg[c], errors="coerce").astype("Int64")
    if asg[["slot_id", "run_order"]].isna().any().any():
        bad = asg[asg.slot_id.isna() | asg.run_order.isna()]
        raise SystemExit(f"ERROR the assignment tab has non-numeric slot_id "
                         f"or run_order:\n{bad.head(10).to_string(index=False)}")
    for t in TARGETS:
        asg[f"well_{t}"] = asg[f"well_{t}"].astype(str).str.strip().str.upper()

    if spec.code.duplicated().any():
        d = spec[spec.code.duplicated(keep=False)]
        raise SystemExit(f"ERROR sample_code is not unique in the specimens "
                         f"tab:\n{d[['sample_code']].to_string(index=False)}")

    asg["is_replicate"] = asg.role.str.upper().ne("PRIMARY")
    asg["replicate_of_plate"] = (asg.role.str.extract(REPLICATE_RE, expand=False)
                                 .str.upper())

    # ---- slot against wells ------------------------------------------------
    expected_wells = asg.slot_id.map(lambda s: slot_wells(int(s)))
    slot_bad = []
    for t, idx in zip(TARGETS, range(3)):
        got = asg[f"well_{t}"]
        want = expected_wells.map(lambda w, i=idx: w[i])
        off = asg.loc[got.ne(want), ["plate_id", "slot_id", f"well_{t}"]].copy()
        if len(off):
            off["expected"] = want[got.ne(want)]
            off["target"] = t
            off = off.rename(columns={f"well_{t}": "recorded"})
            slot_bad.append(off)
    if slot_bad:
        bad = pd.concat(slot_bad, ignore_index=True)
        print(f"ERROR {len(bad)} wells in the assignment tab do not match the "
              f"Layout A mapping their slot implies:", file=sys.stderr)
        print(bad.to_string(index=False), file=sys.stderr)
        print("  slot_id and the well columns are typed independently, so one "
              "of the two is a transcription error. Correct the sample sheet.",
              file=sys.stderr)
        return 4

    # ---- one specimen per well --------------------------------------------
    plate_map = pd.concat(
        [asg[["plate_id", "slot_id", "sample_code", "code", "is_replicate",
              f"well_{t}"]].rename(columns={f"well_{t}": "well"}).assign(target=t)
         for t in TARGETS],
        ignore_index=True).sort_values(["plate_id", "well"])
    dup = plate_map[plate_map.duplicated(["plate_id", "well"], keep=False)]
    if len(dup):
        print(f"ERROR the assignment puts more than one specimen in the same "
              f"well:\n{dup.to_string(index=False)}", file=sys.stderr)
        return 4

    # ---- specimens against assignment, both directions --------------------
    assigned = set(asg.code)
    described = set(spec.code)
    run_undescribed = sorted(assigned - described)
    described_unrun = sorted(described - assigned)
    analysed = assigned & described

    # The status the tracker used to carry by hand, derived instead.
    reused = set(spec.loc[spec.specimen_id.duplicated(keep=False), "code"])
    asg["specimen_status"] = [
        "no specimen record" if c not in described
        else "specimen_id reused by another specimen - keyed on sample_code"
        if c in reused else "OK"
        for c in asg.code]

    n_primary = int((~asg.is_replicate).sum())
    rec = pd.DataFrame([
        {"stage": "slots in the assignment", "n": len(asg), "detail": ""},
        {"stage": "inter-plate replicate slots",
         "n": int(asg.is_replicate.sum()), "detail": ""},
        {"stage": "specimens assigned to a plate", "n": len(assigned),
         "detail": ""},
        {"stage": "specimens described in the specimens tab", "n": len(described),
         "detail": ""},
        {"stage": "assigned and described -- analysed", "n": len(analysed),
         "detail": ""},
        {"stage": "assigned but not described", "n": len(run_undescribed),
         "detail": ",".join(run_undescribed)},
        {"stage": "described but not assigned", "n": len(described_unrun),
         "detail": ",".join(described_unrun)},
        {"stage": "specimen_id reused across specimens", "n": len(reused),
         "detail": ",".join(sorted(reused))},
    ])
    rec.to_csv(args.out_reconciliation, index=False)

    spec.drop(columns="code").to_csv(args.out_specimens, index=False)
    asg.drop(columns="code").to_csv(args.out_assignment, index=False)
    plate_map.to_csv(args.out_plate_map, index=False)

    print(f"sample sheet: {len(asg)} slots on {asg.plate_id.nunique()} plates "
          f"({n_primary} specimens + {len(asg) - n_primary} inter-plate "
          f"replicates)")
    print(f"slot-to-well mapping consistent for all {len(asg)} slots")
    print("\nspecimen reconciliation")
    print(rec[["stage", "n"]].to_string(index=False))

    if described_unrun:
        print(f"\n[!] {len(described_unrun)} specimens are described but were "
              f"never run: {described_unrun[:10]}"
              + (" ..." if len(described_unrun) > 10 else "")
              + "\n    Either a tube was not pulled or a plate in the design "
                "was not run. They cannot contribute to any result.")

    if run_undescribed:
        print(f"\n[!] {len(run_undescribed)} specimens were run but have no "
              f"record in the specimens tab: {run_undescribed[:10]}"
              + (" ..." if len(run_undescribed) > 10 else ""))
        if len(run_undescribed) > args.max_unmatched:
            print(f"\nERROR that is more than --max-unmatched "
                  f"({args.max_unmatched}). A handful of specimens with no "
                  f"record is a documented gap; this many means the two tabs "
                  f"describe different specimens, or sample_code does not "
                  f"line up between them. Analysing the remainder would "
                  f"produce a clean-looking run on an arbitrary subset.",
                  file=sys.stderr)
            return 4
        print(f"    Within --max-unmatched ({args.max_unmatched}). They stay "
              f"in the assignment so the plate-map gate can see their wells, "
              f"and leave the analysis at CALL_TARGETS under "
              f"restrict_to_metadata.")

    print(f"\ndenominator chain: {len(assigned)} run -> {len(run_undescribed)} "
          f"without a record -> {len(analysed)} analysable")
    return 0


if __name__ == "__main__":
    sys.exit(main())
