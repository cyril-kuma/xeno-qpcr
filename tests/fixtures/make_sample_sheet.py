#!/usr/bin/env python3
"""Build the synthetic sample sheet the `test` profile validates against.

Mirrors the real `data/sample_sheet.xlsx`: a `specimens` tab and an
`assignment` tab, each with its header below a one-line note, so
IMPORT_SAMPLE_SHEET is exercised on the same shape.

The fixture deliberately includes one specimen that was run but has no record
in the specimens tab, so CI covers the reconciliation path rather than only
the clean case.

Run after make_plates.py, from the repository root.
"""
from pathlib import Path
import sys

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))
from randomise import slot_wells  # noqa: E402

TARGETS = ("Pf", "Pm", "Po")
REPLICATE_CODE = "TS005"
#: Run but intentionally absent from the specimens tab.
UNDESCRIBED = "TS023"

NOTES = {
    "specimens": ("One row per specimen. Everything known about the mosquito, "
                  "independent of where it was run. Keyed on sample_code."),
    "assignment": ("One row per plate slot. A slot fixes three wells, one per "
                   "target block. Specimen attributes are NOT repeated here; "
                   "they live in the specimens tab."),
}

codes = [f"TS{i:03d}" for i in range(1, 24)]
n = len(codes)
zones = ["Forest", "Northern_Savannah", "Coastal_Savannah"]
sites = ["SiteA", "SiteB", "SiteC"]
methods = ["LTC", "Prokopack_Aspiration", "HLC"]
dates = ["05/05/2022", "18/05/2024", "02/06/2023"]
siblings = ["Anopheles_coluzzii", "Anopheles_gambiae_s.s",
            "Anopheles_arabiensis"]

spec = pd.DataFrame({
    "sample_code": codes,
    "specimen_id": [f"SYN{i:05d}" for i in range(1, n + 1)],
    "sample_id": [f"{c}SYN{i:05d}" for i, c in enumerate(codes, start=1)],
    "species": "Anopheles_gambiae_s.l",
    "sibling_species": [siblings[i % 3] for i in range(n)],
    "feeding_status": "Blood_fed",
    "collection_date": [dates[i % 3] for i in range(n)],
    "collection_time": "06:00am-9:00am",
    "collection_location": [sites[i % 3] for i in range(n)],
    "bioclimatic_zone": [zones[i % 3] for i in range(n)],
    "region": [["Ashanti", "Savannah", "Central"][i % 3] for i in range(n)],
    "latitude": "6.700",
    "longitude": "-1.620",
    "collection_context": "Indoor",
    "collection_method": [methods[i % 3] for i in range(n)],
    "specimen_sex": "Female",
    "specimen_type": "mosquito_engorged_abdomen",
    "specimen_storage": "Ethanol",
    "sample_type": "Wild_adult_caught",
})

# The specimens tab is written to metadata.csv as well, because the plate
# generator and the randomiser both take a plain specimen list.
spec.to_csv(ROOT / "tests" / "data" / "metadata.csv", index=False)

A = codes[:12]
B = codes[12:] + [REPLICATE_CODE]    # 11 specimens + one inter-plate replicate
rows = []
for plate, order, plate_codes in (("R01", 1, A), ("R02", 2, B)):
    for slot, code in enumerate(plate_codes, start=1):
        wells = slot_wells(slot)
        replicate = plate == "R02" and slot == len(B)
        s = spec[spec.sample_code == code]
        rows.append({
            "plate_id": plate, "run_order": order, "slot_id": slot,
            "well_Pf": wells[0], "well_Pm": wells[1], "well_Po": wells[2],
            "sample_code": code,
            "specimen_id": s.specimen_id.iloc[0] if len(s) else None,
            "stratum": (f"{s.bioclimatic_zone.iloc[0]}|"
                        f"{s.collection_location.iloc[0]}|"
                        f"{s.collection_date.iloc[0][-4:]}|"
                        f"{s.collection_method.iloc[0]}") if len(s) else "UNKNOWN",
            "role": ("INTER-PLATE REPLICATE (of plate R01)" if replicate
                     else "PRIMARY"),
            "first_run_plate": "not run in first round",
        })
asg = pd.DataFrame(rows)

# Remove one specimen's record, keeping its slot, so the reconciliation has
# an assigned-but-undescribed specimen to report.
spec_out = spec[spec.sample_code != UNDESCRIBED]

wb = Workbook()
for i, (name, df) in enumerate((("specimens", spec_out), ("assignment", asg))):
    ws = wb.active if i == 0 else wb.create_sheet()
    ws.title = name
    ws.append([NOTES[name]])
    ws.cell(1, 1).font = Font(italic=True, size=9)
    ws.append(list(df.columns))
    for c in range(1, len(df.columns) + 1):
        ws.cell(2, c).font = Font(bold=True)
        ws.cell(2, c).alignment = Alignment(horizontal="center")
    for _, r in df.iterrows():
        ws.append([None if pd.isna(v) else v for v in r.tolist()])
    ws.freeze_panes = "A3"

out = ROOT / "tests" / "data" / "sample_sheet.xlsx"
wb.save(out)
print(f"wrote {out}: {len(spec_out)} specimens, {len(asg)} slots "
      f"({UNDESCRIBED} run but undescribed)")
