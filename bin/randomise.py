#!/usr/bin/env python3
"""Allocate specimens to plates and slots — the design-time first stage.

This is where the study's plate layout came from, and it runs once, before any
plate is set up. Its output is a sample-to-plate assignment: which specimen
occupies which slot on which plate, and therefore which three wells it
occupies in the three target blocks.

WHY RANDOMISE
-------------
Plate-level technical variation is real — a reagent lot, a run, a thermal
block — and if plate is correlated with collection site or bioclimatic zone,
that variation is indistinguishable from a biological difference between
zones. Assay design §8 therefore requires that specimens be distributed so
that plate is independent of site, zone, collection method and year. The
allocation this script produces is checked for exactly that property by
RANDOMISATION_CHECK, and the check is reported with the results.

ALGORITHM
---------
Stratified systematic allocation, which spreads each stratum across plates
more evenly than complete randomisation does at this sample size:

  1. A stratum is zone x site x year x collection method.
  2. Specimens are shuffled within each stratum.
  3. Strata are ordered largest first, ties broken at random.
  4. Specimens are dealt across plates in laps, with a fresh random
     permutation of the plates each lap, so no stratum lands repeatedly on
     the same plate.
  5. Slots left spare after the primary specimens are filled with a second
     aliquot of a specimen already placed on a different plate — the
     inter-plate technical replicates, drawn one each from the largest strata.
  6. Slot order within each plate is randomised, so position on the plate
     carries no information either.
  7. Each slot maps to its three wells by the Layout A block structure
     (assay design §6).

THE FROZEN ALLOCATION
---------------------
The allocation used in this study was produced at seed 20261001 and is
recorded in the `assignment` tab of `data/sample_sheet.xlsx`, which is what
the plates were set up from and what IMPORT_SAMPLE_SHEET reads.
Re-running this script does not alter the
analysis, and it will not reproduce the recorded allocation bit-for-bit
unless the inputs and the NumPy generator match those used at the time. To
allocate a new set of specimens, run this script, commit its output, and
write it into a sample sheet before setting up any plate.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

TARGETS = ("Pf", "Pm", "Po")
#: Layout A sample slots: slots 1-4 occupy rows E-H of the block's third
#: column, slots 5-12 occupy rows A-H of its fourth (assay design §6).
ROWS = "ABCDEFGH"


def read_specimens_tab(book: Path) -> pd.DataFrame:
    """The specimens tab of a sample sheet, whose header sits below a note."""
    probe = pd.read_excel(book, "specimens", engine="calamine",
                          header=None, nrows=4)
    header = 0
    for i in range(min(4, len(probe))):
        if "sample_code" in [str(v).strip() for v in probe.iloc[i].tolist()]:
            header = i
            break
    return pd.read_excel(book, "specimens", engine="calamine", dtype=str,
                         header=header)


def slot_wells(slot: int) -> tuple[str, str, str]:
    """The three wells a slot occupies, one per target block."""
    if slot <= 4:
        row, col = ROWS[slot + 3], 3
    else:
        row, col = ROWS[slot - 5], 4
    return f"{row}{col}", f"{row}{col + 4}", f"{row}{col + 8}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--specimens", required=True, type=Path,
                    help="the specimens to allocate: a CSV, or a sample "
                         "sheet whose 'specimens' tab holds them")
    ap.add_argument("--first-round-slots", type=Path,
                    help="optional first-round plate layout with 'Sample "
                         "Name' and 'plate' columns. Supplies first_run_plate "
                         "and carries over specimens that were run but are "
                         "absent from the metadata")
    ap.add_argument("--n-plates", type=int, default=60)
    ap.add_argument("--slots-per-plate", type=int, default=12,
                    help="Layout A carries 12 specimens per plate (assay "
                         "design §6)")
    ap.add_argument("--seed", type=int, required=True,
                    help="record this in the laboratory notebook alongside "
                         "the date of randomisation")
    ap.add_argument("--out-assignment", required=True, type=Path)
    ap.add_argument("--out-balance", required=True, type=Path)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    if args.specimens.suffix.lower() in (".xlsx", ".xlsm"):
        meta = read_specimens_tab(args.specimens)
    else:
        meta = pd.read_csv(args.specimens, dtype=str)
    for col in ("sample_code", "specimen_id"):
        if col not in meta.columns:
            print(f"ERROR {args.specimens.name} has no '{col}' column; "
                  f"columns: {list(meta.columns)}", file=sys.stderr)
            return 4

    first: dict[str, str] = {}
    extra: list[str] = []
    if args.first_round_slots and args.first_round_slots.exists():
        fr = pd.read_csv(args.first_round_slots, dtype=str)
        fr = fr.dropna(subset=["Sample Name"])
        fr["sample_code"] = fr["Sample Name"].str.replace("_", "", regex=False)
        first = fr.set_index("sample_code")["plate"].to_dict()
        extra = sorted(set(fr.sample_code) - set(meta.sample_code))

    S = pd.concat(
        [meta, pd.DataFrame({"sample_code": extra, "sample_id": extra})],
        ignore_index=True)
    S["year"] = S.collection_date.str[-4:]
    S["stratum"] = (S.bioclimatic_zone.fillna("NO_METADATA") + "|"
                    + S.collection_location.fillna("NA") + "|"
                    + S.year.fillna("NA") + "|"
                    + S.collection_method.fillna("NA"))
    S["metadata_status"] = np.where(
        S.bioclimatic_zone.isna(), "MISSING - resolve before analysis",
        np.where(S.specimen_id.duplicated(keep=False),
                 "specimen_id reused by another specimen - key on sample_code",
                 "OK"))
    S["first_run_plate"] = S.sample_code.map(first).fillna(
        "not run in first round")

    capacity = args.n_plates * args.slots_per_plate
    if len(S) > capacity:
        print(f"ERROR {len(S)} specimens exceed the capacity of "
              f"{args.n_plates} plates x {args.slots_per_plate} slots "
              f"({capacity}). Increase --n-plates to "
              f"{-(-len(S) // args.slots_per_plate)}.", file=sys.stderr)
        return 4

    # ---- steps 2-4: shuffle within strata, order strata, deal in laps ----
    sizes = S.stratum.value_counts()
    strata = sorted(sizes.index, key=lambda s: (-sizes[s], rng.random()))
    order: list[int] = []
    for st in strata:
        idx = S.index[S.stratum == st].to_numpy().copy()
        rng.shuffle(idx)
        order += list(idx)

    plate = np.empty(len(order), int)
    for lap in range(int(np.ceil(len(order) / args.n_plates))):
        perm = rng.permutation(args.n_plates)
        lo = lap * args.n_plates
        for k in range(len(order[lo:lo + args.n_plates])):
            plate[lo + k] = perm[k]
    S.loc[order, "plate_idx"] = plate
    S["role"] = "PRIMARY"

    # ---- step 5: spare slots become inter-plate technical replicates ----
    cnt = S.plate_idx.value_counts()
    short = list(cnt[cnt < args.slots_per_plate].index)
    rng.shuffle(short)
    reps = []
    for st, pl in zip(strata[:len(short)], short):
        cand = S[(S.stratum == st) & (S.plate_idx != pl) & (S.role == "PRIMARY")]
        if cand.empty:
            continue
        r = cand.loc[rng.choice(cand.index)].copy()
        origin = int(S.loc[S.sample_code == r.sample_code, "plate_idx"].iloc[0])
        r["plate_idx"] = pl
        r["role"] = f"INTER-PLATE REPLICATE (of plate R{origin + 1:02d})"
        reps.append(r)
    if reps:
        S = pd.concat([S, pd.DataFrame(reps)], ignore_index=True)

    # ---- steps 6-7: randomise slot within plate, then map slot to wells ----
    S["slot_id"] = 0
    for _pl, d in S.groupby("plate_idx"):
        S.loc[d.index, "slot_id"] = rng.permutation(len(d)) + 1
    S["plate_id"] = S.plate_idx.astype(int).map(lambda x: f"R{x + 1:02d}")
    S["run_order"] = S.plate_idx.astype(int) + 1
    S[list(f"well_{t}" for t in TARGETS)] = S.slot_id.apply(
        lambda s: pd.Series(slot_wells(int(s))))
    S["is_replicate"] = S.role.ne("PRIMARY")
    S["replicate_of_plate"] = S.role.str.extract(
        r"INTER-PLATE REPLICATE \(of plate ([A-Za-z0-9-]+)\)", expand=False)

    cols = ["plate_id", "run_order", "slot_id", "well_Pf", "well_Pm", "well_Po",
            "specimen_id", "sample_code", "sample_id", "collection_location",
            "bioclimatic_zone", "region", "sibling_species", "collection_date",
            "collection_context", "collection_method", "specimen_storage",
            "year", "stratum", "role", "is_replicate", "replicate_of_plate",
            "first_run_plate", "metadata_status"]
    S = (S.sort_values(["run_order", "slot_id"]).reset_index(drop=True)
          [[c for c in cols if c in S.columns]])
    S.to_csv(args.out_assignment, index=False)

    # ---- the balance the randomisation exists to achieve ----
    checked = ("bioclimatic_zone", "collection_method", "collection_location",
               "year", "stratum", "sibling_species")
    absent = [c for c in checked if c not in S.columns]
    if absent:
        print(f"WARNING cannot check balance against {absent}: the metadata "
              f"has no such column, so the allocation is unverified on it",
              file=sys.stderr)
    rows = []
    for c in checked:
        if c not in S.columns:
            continue
        t = pd.crosstab(S.plate_id, S[c].fillna("NA"))
        chi2, p, _dof, _exp = chi2_contingency(t)
        rows.append({"variable": c, "chi2": round(float(chi2), 2),
                     "p_value": round(float(p), 4),
                     "max_level_share": round(
                         float(t.div(t.sum(axis=1), axis=0).max().max()), 3)})
    if "first_run_plate" in S.columns:
        t = pd.crosstab(S.plate_id, S.first_run_plate).drop(
            columns="not run in first round", errors="ignore")
        if t.shape[1]:
            rows.append({"variable": "first_run_plate", "chi2": np.nan,
                         "p_value": np.nan,
                         "max_level_share": int(t.max().max())})
    bal = pd.DataFrame(rows)
    bal.to_csv(args.out_balance, index=False)

    n_rep = int(S.is_replicate.sum())
    print(f"seed {args.seed}: {len(S)} slots on {S.plate_id.nunique()} plates "
          f"({len(S) - n_rep} specimens + {n_rep} inter-plate replicates)")
    print(f"slots per plate: {S.groupby('plate_id').size().value_counts().to_dict()}")
    print("\nallocation balance")
    print(bal.to_string(index=False))
    print("\nRecord the seed, the date and this script's version in the "
          "laboratory notebook, and commit the assignment before setting up "
          "any plate.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
