#!/usr/bin/env python3
"""Cross-reactivity verdict from the validation plates.

The co-detection exclusion rule in call_targets.py treats a Po positive
co-occurring with Pf as the Po reaction amplifying Pf template. That rule
rests on an empirical claim about the assay's specificity, and the place to
establish it is a validation plate: known-composition material run through all
three reactions, where a signal in the wrong channel has no other explanation.

This stage produces the verdict rather than assuming it. Three states:

  clean            every target assay was negative on material not containing
                   that species, and the uninfected-host controls were negative
  cross_reactive   at least one target amplified material not containing that
                   species -- the co-detection rule is supported by direct
                   evidence, and the detail names the channels involved
  not_assessed     no validation plate was supplied

`not_assessed` is not a failure. It records that specificity is supported by
the per-plate uninfected-host negative controls and by the published primer
design rather than by a validation plate in this run, and call_targets.py
carries that state into the output so the report states it instead of leaving
the reader to assume a validation was done.

Validation plates are identified by `plate_class == "validation"`, which
ingest derives from the plate identifier.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from well_roles import NEGATIVE_HOST, POSITIVE

ap = argparse.ArgumentParser()
ap.add_argument("--wells", required=True, type=Path)
ap.add_argument("--targets", default="Pf,Pm,Po")
ap.add_argument("--cq-positive-max", type=float, default=38.0)
ap.add_argument("--out-verdict", required=True, type=Path)
ap.add_argument("--out-detail", required=True, type=Path)
a = ap.parse_args()

targets = [t.strip() for t in a.targets.split(",") if t.strip()]
w = pd.read_csv(a.wells)
val = w[w.plate_class.eq("validation")].copy()

rows: list[dict] = []

if val.empty:
    verdict = {
        "verdict": "not_assessed",
        "n_validation_plates": 0,
        "detail": ("no validation plate was supplied, so the absence of "
                   "cross-reactivity with the uninfected Anopheles genome "
                   "rests on the per-plate uninfected-host negative controls "
                   "alone, together with the published primer design"),
    }
    rows.append({"check": "validation_plates", "n": 0,
                 "detail": "none supplied"})
else:
    amp = val.cq.notna() & (val.cq <= a.cq_positive_max)
    val["amplified"] = amp

    # Uninfected-host gDNA amplifying in any Plasmodium channel is direct
    # evidence of cross-reactivity with the mosquito genome.
    host = val[val.role.eq(NEGATIVE_HOST)]
    host_hits = host[host.amplified]

    # A positive control for one species amplifying in another species'
    # channel is direct evidence of cross-reactivity between target assays.
    pos = val[val.role.eq(POSITIVE) & val.control_target.notna()]
    off_target = pos[pos.amplified & pos.control_target.ne(pos.target)]

    for t in targets:
        rows.append({"check": f"host_gdna_in_{t}_channel",
                     "n": int(host_hits.target.eq(t).sum()),
                     "detail": "uninfected Anopheles gDNA amplifying as " + t})
    for _, r in off_target.iterrows():
        rows.append({"check": "off_target_positive_control", "n": 1,
                     "detail": (f"{r.sample_name} (target {r.control_target}) "
                                f"amplified in the {r.target} channel at "
                                f"Cq {r.cq} on plate {r.plate}")})

    n_bad = len(host_hits) + len(off_target)
    plates = sorted(val.plate.unique())
    if n_bad:
        chans = sorted(set(host_hits.target) | set(off_target.target))
        verdict = {
            "verdict": "cross_reactive",
            "n_validation_plates": len(plates),
            "detail": (f"{n_bad} wells on validation plates "
                       f"{', '.join(plates)} amplified material not "
                       f"containing the target species; channels involved: "
                       f"{', '.join(chans)}. The co-detection exclusion rule "
                       f"is supported by direct evidence."),
        }
    else:
        verdict = {
            "verdict": "clean",
            "n_validation_plates": len(plates),
            "detail": (f"no cross-reactivity on validation plates "
                       f"{', '.join(plates)}: uninfected-host gDNA was "
                       f"negative in all three channels and every positive "
                       f"control amplified only in its own channel. Any "
                       f"co-detection in the field specimens therefore needs "
                       f"an explanation other than channel cross-talk."),
        }
    rows.append({"check": "validation_plates", "n": len(plates),
                 "detail": ",".join(plates)})

a.out_verdict.write_text(json.dumps(verdict, indent=2))
pd.DataFrame(rows).to_csv(a.out_detail, index=False)
print(f"cross-reactivity: {verdict['verdict']}")
print(f"  {verdict['detail']}")
