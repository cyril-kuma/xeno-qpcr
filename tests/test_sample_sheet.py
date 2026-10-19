"""The sample sheet: its two checks, and the reconciliation it reports."""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest
from openpyxl import Workbook

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

from randomise import slot_wells  # noqa: E402

TARGETS = ("Pf", "Pm", "Po")

SPEC_COLS = ["sample_code", "specimen_id", "bioclimatic_zone",
             "collection_location", "collection_method", "collection_date"]
ASSIGN_COLS = ["plate_id", "run_order", "slot_id", "well_Pf", "well_Pm",
               "well_Po", "sample_code", "specimen_id", "stratum", "role",
               "first_run_plate"]


def _sheet(path, specimens, assignment, note=True):
    """Write a two-tab sample sheet, header below a one-line note."""
    wb = Workbook()
    for i, (name, df) in enumerate((("specimens", specimens),
                                    ("assignment", assignment))):
        ws = wb.active if i == 0 else wb.create_sheet()
        ws.title = name
        if note:
            ws.append([f"note about the {name} tab"])
        ws.append(list(df.columns))
        for _, r in df.iterrows():
            ws.append([None if pd.isna(v) else v for v in r.tolist()])
    wb.save(path)


def _cycle(values, n):
    return [values[i % len(values)] for i in range(n)]


def _data(n_specimens=6, n_slots=6, drop_records=(), well_override=None,
          duplicate_specimen_id=False):
    codes = [f"TS{i:03d}" for i in range(1, max(n_specimens, n_slots) + 1)]
    ids = [f"SYN{i}" for i in range(1, n_specimens + 1)]
    if duplicate_specimen_id and n_specimens >= 2:
        ids[1] = ids[0]
    spec = pd.DataFrame({
        "sample_code": codes[:n_specimens],
        "specimen_id": ids,
        "bioclimatic_zone": _cycle(["Forest", "Northern_Savannah"], n_specimens),
        "collection_location": _cycle(["SiteA", "SiteB"], n_specimens),
        "collection_method": _cycle(["LTC", "HLC"], n_specimens),
        "collection_date": _cycle(["05/05/2022", "18/05/2024"], n_specimens),
    })
    spec = spec[~spec.sample_code.isin(drop_records)]

    rows = []
    for slot in range(1, n_slots + 1):
        w = list(slot_wells(slot))
        if well_override and slot == well_override[0]:
            w[well_override[1]] = well_override[2]
        rows.append({"plate_id": "R01", "run_order": 1, "slot_id": slot,
                     "well_Pf": w[0], "well_Pm": w[1], "well_Po": w[2],
                     "sample_code": codes[slot - 1],
                     "specimen_id": f"SYN{slot}", "stratum": "s", "role": "PRIMARY",
                     "first_run_plate": "not run in first round"})
    return spec, pd.DataFrame(rows)


def _run(specimens, assignment, extra=(), note=True):
    with tempfile.TemporaryDirectory() as d:
        sheet = os.path.join(d, "sample_sheet.xlsx")
        _sheet(sheet, specimens, assignment, note=note)
        out = {k: os.path.join(d, f"{k}.csv")
               for k in ("specimens", "assignment", "plate_map", "reconciliation")}
        r = subprocess.run(
            [sys.executable, str(BIN / "import_sample_sheet.py"),
             "--sample-sheet", sheet,
             "--out-specimens", out["specimens"],
             "--out-assignment", out["assignment"],
             "--out-plate-map", out["plate_map"],
             "--out-reconciliation", out["reconciliation"], *extra],
            capture_output=True, text=True,
            env={**os.environ, "PYTHONPATH": str(BIN)})
        loaded = {k: (pd.read_csv(v) if os.path.exists(v) else None)
                  for k, v in out.items()}
        return r, loaded


def test_a_clean_sheet_imports():
    r, out = _run(*_data())
    assert r.returncode == 0, r.stderr
    rec = dict(zip(out["reconciliation"].stage, out["reconciliation"].n))
    assert rec["specimens assigned to a plate"] == 6
    assert rec["assigned and described -- analysed"] == 6
    assert rec["assigned but not described"] == 0
    assert rec["described but not assigned"] == 0
    # three wells per slot
    assert len(out["plate_map"]) == 18
    assert set(out["plate_map"].target) == set(TARGETS)


def test_header_found_with_or_without_the_note_row():
    for note in (True, False):
        r, out = _run(*_data(), note=note)
        assert r.returncode == 0, f"note={note}: {r.stderr}"
        assert len(out["assignment"]) == 6


def test_a_mistyped_well_is_caught_against_the_slot():
    """slot_id and the well columns are typed independently but related by the
    Layout A mapping, so one contradicts the other. This is the check that
    replaced comparing two generated tabs, which would have been circular."""
    spec, asg = _data(well_override=(1, 1, "E8"))     # slot 1 Pm should be E7
    r, _out = _run(spec, asg)
    assert r.returncode == 4
    assert "do not match the Layout A mapping" in r.stderr
    assert "E8" in r.stderr and "E7" in r.stderr


def test_two_specimens_in_one_well_is_caught():
    """A duplicated slot: both rows agree with the Layout A mapping, so the
    slot check passes and the two specimens collide in the same wells."""
    spec, asg = _data()
    asg.loc[1, "slot_id"] = 1
    asg.loc[1, ["well_Pf", "well_Pm", "well_Po"]] = list(slot_wells(1))
    r, _out = _run(spec, asg)
    assert r.returncode == 4
    assert "more than one specimen in the same well" in r.stderr


def test_an_undescribed_specimen_is_reported_not_dropped():
    """It is physically in a well on a real plate, so the plate-map gate must
    still see it. It leaves the analysis later, at CALL_TARGETS."""
    r, out = _run(*_data(drop_records=("TS003",)))
    assert r.returncode == 0, r.stderr
    rec = out["reconciliation"].set_index("stage")
    assert rec.loc["assigned but not described", "n"] == 1
    assert rec.loc["assigned but not described", "detail"] == "TS003"
    assert rec.loc["assigned and described -- analysed", "n"] == 5
    # still in the assignment and still in the plate map
    assert "TS003" in set(out["assignment"].sample_code)
    assert (out["plate_map"].sample_code == "TS003").sum() == 3
    # and flagged, so the status is derived rather than hand-maintained
    status = dict(zip(out["assignment"].sample_code,
                      out["assignment"].specimen_status))
    assert status["TS003"] == "no specimen record"
    assert status["TS001"] == "OK"


def test_too_many_undescribed_specimens_aborts():
    """A handful of gaps and a broken join look identical to a silent
    intersection: both yield a clean run on whatever survived."""
    spec, asg = _data(n_specimens=6, n_slots=6,
                      drop_records=("TS002", "TS003", "TS004", "TS005"))
    r, _out = _run(spec, asg, extra=("--max-unmatched", "2"))
    assert r.returncode == 4
    assert "more than --max-unmatched" in r.stderr

    # the same sheet passes when the threshold admits it, and still records it
    r, out = _run(spec, asg, extra=("--max-unmatched", "99"))
    assert r.returncode == 0, r.stderr
    assert out["reconciliation"].set_index("stage").loc[
        "assigned but not described", "n"] == 4


def test_a_described_specimen_that_was_never_run_is_reported():
    """The other direction: a tube not pulled, or a plate not run."""
    r, out = _run(*_data(n_specimens=8, n_slots=6))
    assert r.returncode == 0, r.stderr
    rec = out["reconciliation"].set_index("stage")
    assert rec.loc["described but not assigned", "n"] == 2
    assert "never run" in r.stdout


def test_reused_specimen_id_is_derived_not_declared():
    r, out = _run(*_data(duplicate_specimen_id=True))
    assert r.returncode == 0, r.stderr
    rec = out["reconciliation"].set_index("stage")
    assert rec.loc["specimen_id reused across specimens", "n"] == 2
    status = dict(zip(out["assignment"].sample_code,
                      out["assignment"].specimen_status))
    assert "keyed on sample_code" in status["TS001"]


def test_duplicate_sample_code_in_specimens_aborts():
    """sample_code is the join key everywhere, so it must be unique."""
    spec, asg = _data()
    spec = pd.concat([spec, spec.iloc[[0]]], ignore_index=True)
    r, _out = _run(spec, asg)
    assert r.returncode != 0
    assert "not unique" in r.stderr


def test_replicate_role_is_parsed_into_its_origin_plate():
    spec, asg = _data()
    asg.loc[5, "role"] = "INTER-PLATE REPLICATE (of plate R39)"
    asg.loc[5, "sample_code"] = "TS001"
    asg.loc[5, ["well_Pf", "well_Pm", "well_Po"]] = list(slot_wells(6))
    r, out = _run(spec, asg)
    assert r.returncode == 0, r.stderr
    row = out["assignment"].iloc[5]
    assert bool(row.is_replicate)
    assert row.replicate_of_plate == "R39"
    assert out["reconciliation"].set_index("stage").loc[
        "inter-plate replicate slots", "n"] == 1


@pytest.mark.parametrize("tab,col", [("specimens", "sample_code"),
                                     ("assignment", "plate_id")])
def test_a_missing_required_column_names_what_was_found(tab, col):
    spec, asg = _data()
    if tab == "specimens":
        spec = spec.drop(columns=["bioclimatic_zone"])
    else:
        asg = asg.drop(columns=["stratum"])
    r, _out = _run(spec, asg)
    assert r.returncode != 0
    assert "is missing" in r.stderr and "Columns found" in r.stderr


def test_the_real_sample_sheet_imports_cleanly():
    """The shipped sheet must satisfy every check in this file."""
    sheet = ROOT / "data" / "sample_sheet.xlsx"
    if not sheet.exists():
        pytest.skip("the specimen data is not present in this checkout")
    with tempfile.TemporaryDirectory() as d:
        out = {k: os.path.join(d, f"{k}.csv")
               for k in ("specimens", "assignment", "plate_map", "reconciliation")}
        r = subprocess.run(
            [sys.executable, str(BIN / "import_sample_sheet.py"),
             "--sample-sheet", str(sheet),
             "--out-specimens", out["specimens"],
             "--out-assignment", out["assignment"],
             "--out-plate-map", out["plate_map"],
             "--out-reconciliation", out["reconciliation"]],
            capture_output=True, text=True,
            env={**os.environ, "PYTHONPATH": str(BIN)})
        assert r.returncode == 0, r.stderr
        rec = pd.read_csv(out["reconciliation"]).set_index("stage")
        assert rec.loc["slots in the assignment", "n"] == 720
        assert rec.loc["inter-plate replicate slots", "n"] == 6
        assert rec.loc["specimens assigned to a plate", "n"] == 714
        assert rec.loc["assigned and described -- analysed", "n"] == 701
        assert rec.loc["assigned but not described", "n"] == 13
        assert rec.loc["described but not assigned", "n"] == 0
        assert len(pd.read_csv(out["plate_map"])) == 2160
