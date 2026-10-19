"""The gate: the plates that ran must be the plates the design specified."""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

from randomise import slot_wells  # noqa: E402

TARGETS = ("Pf", "Pm", "Po")


def _fixture(swap=False, drop=False):
    """Two plates of three specimens each, and the assignment describing them.

    `swap` exchanges two specimens between wells on one plate: the calls stay
    internally consistent and attach to the wrong specimen, which is the
    failure no later stage can detect. `drop` removes a well the assignment
    expects.
    """
    rows, exp = [], []
    for plate, codes in (("R01", ["S1", "S2", "S3"]),
                         ("R02", ["S4", "S5", "S6"])):
        for slot, code in enumerate(codes, start=1):
            wells = slot_wells(slot)
            for t, well in zip(TARGETS, wells):
                exp.append({"plate_id": plate, "slot_id": slot,
                            "sample_code": code, "code": code,
                            "is_replicate": False, "well": well,
                            "target": t})
                rows.append({"plate": plate, "plate_class": "routine",
                             "role": "sample", "well": well, "target": t,
                             "sample_name": code, "sample_code": code,
                             "cq": 25.0, "cq_conf": 0.9, "amp_status": "Amp"})
    wells_df = pd.DataFrame(rows)
    if swap:
        a = (wells_df.plate == "R01") & (wells_df.sample_code == "S1")
        b = (wells_df.plate == "R01") & (wells_df.sample_code == "S2")
        wells_df.loc[a, ["sample_code", "sample_name"]] = "S2"
        wells_df.loc[b, ["sample_code", "sample_name"]] = "S1"
    if drop:
        wells_df = wells_df[~((wells_df.plate == "R02")
                              & (wells_df.target == "Po"))]
    return wells_df, pd.DataFrame(exp)


def _run(wells_df, exp_df, extra=()):
    with tempfile.TemporaryDirectory() as d:
        w = os.path.join(d, "wells_merged.csv")
        e = os.path.join(d, "plate_map.csv")
        ow = os.path.join(d, "wells.csv")
        orp = os.path.join(d, "report.csv")
        wells_df.to_csv(w, index=False)
        exp_df.to_csv(e, index=False)
        r = subprocess.run(
            [sys.executable, str(BIN / "validate_plate_maps.py"),
             "--wells", w, "--plate-map", e, "--out-wells", ow,
             "--out-report", orp, *extra],
            capture_output=True, text=True)
        rep = pd.read_csv(orp) if os.path.exists(orp) else None
        return r, rep


def test_a_matching_plate_set_passes():
    r, rep = _run(*_fixture())
    assert r.returncode == 0, r.stderr
    got = dict(zip(rep.check, rep.n))
    assert got["wells_matched"] == 18
    assert got["wells_misplaced"] == 0
    assert got["wells_unexpected"] == 0
    assert got["wells_absent"] == 0


def test_two_specimens_swapped_between_wells_aborts():
    """The failure the gate exists for. Every call is valid; each is simply
    attached to the wrong specimen, and therefore to the wrong zone, site
    and collection method."""
    r, rep = _run(*_fixture(swap=True))
    assert r.returncode == 6
    assert dict(zip(rep.check, rep.n))["wells_misplaced"] == 6
    assert "do not match the assignment" in r.stderr


def test_a_swap_within_one_plate_is_named_as_such():
    r, rep = _run(*_fixture(swap=True))
    assert dict(zip(rep.check, rep.n))["wells_wrong_slot_same_plate"] == 6
    assert "wrong slot" in r.stderr


def test_a_missing_well_aborts():
    r, rep = _run(*_fixture(drop=True))
    assert r.returncode == 6
    assert dict(zip(rep.check, rep.n))["wells_absent"] == 3


def test_the_override_records_the_mismatch_and_continues():
    r, rep = _run(*_fixture(swap=True),
                  extra=("--allow-plate-map-mismatch",))
    assert r.returncode == 0
    assert dict(zip(rep.check, rep.n))["wells_misplaced"] == 6
    assert "WARNING" in r.stderr


def test_validation_plate_wells_are_exempt():
    """Validation plates carry known-composition material, which the
    specimen assignment does not describe."""
    wells_df, exp_df = _fixture()
    extra = pd.DataFrame([{
        "plate": "V01", "plate_class": "validation", "role": "sample",
        "well": "H4", "target": "Pf", "sample_name": "VAL_1",
        "sample_code": "VAL1", "cq": 20.0, "cq_conf": 0.9,
        "amp_status": "Amp"}])
    r, rep = _run(pd.concat([wells_df, extra], ignore_index=True), exp_df)
    assert r.returncode == 0, r.stderr
    assert dict(zip(rep.check, rep.n))["wells_unexpected"] == 0


def test_slot_wells_match_the_layout_a_block_structure():
    """Assay design §6: slots 1-4 in rows E-H of column 3, slots 5-12 in rows
    A-H of column 4, with the three target blocks four columns apart. This is
    the mapping the tracker's Sample_Assignment records."""
    assert slot_wells(1) == ("E3", "E7", "E11")
    assert slot_wells(4) == ("H3", "H7", "H11")
    assert slot_wells(5) == ("A4", "A8", "A12")
    assert slot_wells(12) == ("H4", "H8", "H12")
    seen = [w for s in range(1, 13) for w in slot_wells(s)]
    assert len(set(seen)) == 36
