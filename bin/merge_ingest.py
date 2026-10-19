#!/usr/bin/env python3
"""Concatenate per-plate well tables and flatten the ingest audits."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--wells", nargs="+", required=True)
ap.add_argument("--audits", nargs="+", required=True)
ap.add_argument("--out-wells", required=True, type=Path)
ap.add_argument("--out-audit", required=True, type=Path)
a = ap.parse_args()

w = pd.concat([pd.read_csv(f) for f in sorted(a.wells)], ignore_index=True)
w.to_csv(a.out_wells, index=False)

rows = []
for f in sorted(a.audits):
    j = json.loads(Path(f).read_text())
    rows.append({"plate": j["plate"], "plate_class": j["plate_class"],
                 "file": j["file"], "header_row": j["header_row"],
                 "n_wells": j["n_wells"], "targets": ",".join(j["targets"]),
                 "status": j["status"],
                 "n_drifted_cells": len(j["ingest_drift"]),
                 "drift_detail": json.dumps(j["ingest_drift"]) if j["ingest_drift"] else ""})
aud = pd.DataFrame(rows)
aud.to_csv(a.out_audit, index=False)

by_class = aud.plate_class.value_counts().to_dict()
print(f"merged {len(w)} wells from {len(aud)} plates {by_class}")

dup = aud.plate[aud.plate.duplicated()].tolist()
if dup:
    print(f"ERROR plate identifiers are not unique: {sorted(set(dup))}. Two "
          f"exports resolved to the same plate, so their wells would be "
          f"pooled. Give a repeat run a suffix, e.g. R07-2.", file=sys.stderr)
    sys.exit(3)

bad = aud[aud.status != "ok"]
if len(bad):
    print(f"ERROR {len(bad)} plates with ingest drift", file=sys.stderr)
    print(bad[["plate", "n_drifted_cells"]].to_string(index=False), file=sys.stderr)
    sys.exit(2)
