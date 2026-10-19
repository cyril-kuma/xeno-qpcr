"""Regression tests locking the calling rule and the faults it guards against."""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

from call_targets import (apply_codetection_exclusion, apply_guardband,  # noqa: E402
                          apply_ipc, base_call, normalise_key)

P = {"cq_positive_max": 38.0, "cq_trace_max": 40.0,
     "cqconf_min": 0.50, "cqconf_targets": {"Pm"}}


@pytest.mark.parametrize("amp,cq,tgt,conf,want", [
    ("Inconclusive", 25.0, "Pf", 0.9, "invalid"),
    ("No Amp",       None, "Pf", 0.0, "negative"),
    ("Amp",          29.0, "Pf", 0.9, "positive"),
    ("Amp",          39.0, "Pf", 0.9, "trace"),
    ("Amp",          41.0, "Pf", 0.9, "negative"),
    ("Amp",          38.0, "Pf", 0.9, "positive"),   # boundary is inclusive
    ("Amp",          30.0, "Pm", 0.9, "positive"),
    ("Amp",          30.0, "Pm", 0.2, "trace"),      # Cq-confidence floor
    ("Amp",          30.0, "Pm", None, "trace"),
    ("Amp",          30.0, "Po", 0.1, "positive"),   # floor applies to Pm only
])
def test_base_call(amp, cq, tgt, conf, want):
    assert base_call(amp, cq, tgt, conf, P) == want


def test_no_amp_wins_over_a_present_cq():
    """Some exports carry a Cq on a well the instrument called `No Amp` --
    plate R39 well A4 reads Cq 22.71 with Amp Status `No Amp`. The
    amplification status is authoritative; a Cq alone is not a detection."""
    assert base_call("No Amp", 22.709791, "Pf", 0.139432, P) == "negative"


def test_guardband_reaches_every_non_dominant_target():
    """An earlier implementation demoted Pm but never Po, so Po returned an
    identical count at every guardband. Both must respond."""
    wide = pd.DataFrame({
        "Pf_call": ["positive"] * 4, "Pf_cq": [25.0] * 4,
        "Pm_call": ["positive"] * 4, "Pm_cq": [24.0, 26.0, 29.0, 31.0],
        "Po_call": ["positive"] * 4, "Po_cq": [24.0, 26.0, 29.0, 31.0],
    })
    for t in ("Pm", "Po"):
        counts = {
            gb: (apply_guardband(wide, ["Pf", "Pm", "Po"], "Pf", gb)[f"{t}_call_gb"]
                 == "positive").sum()
            for gb in (0.0, 2.0, 5.0)
        }
        assert len(set(counts.values())) > 1, f"guardband is a no-op for {t}: {counts}"


def test_guardband_never_touches_dominant():
    wide = pd.DataFrame({"Pf_call": ["positive"], "Pf_cq": [25.0],
                         "Pm_call": ["positive"], "Pm_cq": [30.0],
                         "Po_call": ["negative"], "Po_cq": [None]})
    out = apply_guardband(wide, ["Pf", "Pm", "Po"], "Pf", 0.0)
    assert out["Pf_call_gb"].iloc[0] == "positive"


def test_codetection_exclusion_is_stricter_than_guardband_zero():
    """A co-detection whose Cq falls BELOW the dominant target's survives a
    guardband of 0 but must not survive co-detection exclusion."""
    wide = pd.DataFrame({
        "Pf_call": ["positive", "negative"], "Pf_cq": [30.0, None],
        "Pm_call": ["positive", "positive"], "Pm_cq": [28.0, 28.0],
        "Po_call": ["positive", "positive"], "Po_cq": [28.0, 28.0],
    })
    gb = apply_guardband(wide, ["Pf", "Pm", "Po"], "Pf", 0.0)
    assert (gb["Po_call_gb"] == "positive").sum() == 2      # guardband keeps both

    ex = apply_codetection_exclusion(wide, ["Pf", "Pm", "Po"], "Pf", {"Po"})
    assert ex["Po_call"].tolist() == ["excluded_codetection", "positive"]
    assert ex["Pm_call"].tolist() == ["positive", "positive"]   # Pm not listed
    assert ex["Pf_call"].tolist() == ["positive", "negative"]   # dominant untouched


def test_ipc_failure_becomes_indeterminate_never_negative():
    """Without a passing IPC a well with no Cq has not demonstrated that the
    extraction worked, so it cannot be reported as a true negative."""
    wide = pd.DataFrame({
        "ipc_cq": [25.0, 35.0, None],
        "Pf_call": ["negative", "negative", "negative"],
        "Pm_call": ["trace", "trace", "trace"],
        "Po_call": ["positive", "positive", "positive"],
    })
    out = apply_ipc(wide, ["Pf", "Pm", "Po"], 32.0)
    assert out["ipc_pass"].tolist() == [True, False, False]
    # the passing specimen is untouched
    assert out.loc[0, ["Pf_call", "Pm_call", "Po_call"]].tolist() == [
        "negative", "trace", "positive"]
    # the failing ones lose their negative and trace calls
    assert out.loc[1, "Pf_call"] == "indeterminate"
    assert out.loc[2, "Pm_call"] == "indeterminate"
    # a target that amplified did so regardless of what the IPC shows
    assert out.loc[1, "Po_call"] == "positive"
    assert out.loc[2, "Po_call"] == "positive"


def test_indeterminate_is_excluded_from_the_denominator():
    """An IPC failure must leave the denominator, not sit in it as a
    negative -- otherwise the IPC would deflate prevalence rather than
    qualify it."""
    wells = pd.DataFrame({
        "role": ["sample"] * 8,
        "plate_class": ["routine"] * 8,
        "plate": ["R01"] * 8,
        "well": ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4"],
        "sample_code": ["S1"] * 4 + ["S2"] * 4,
        "target": ["Pf", "Pm", "Po", "IPC"] * 2,
        "cq": [25.0, None, None, 24.0, None, None, None, None],
        "cq_conf": [0.9, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.0],
        "amp_status": ["Amp", "No Amp", "No Amp", "Amp"] + ["No Amp"] * 4,
    })
    with tempfile.TemporaryDirectory() as d:
        w = os.path.join(d, "wells.csv")
        wells.to_csv(w, index=False)
        out = [os.path.join(d, f) for f in ("c.csv", "m.csv", "s.csv", "i.csv")]
        r = subprocess.run(
            [sys.executable, str(BIN / "call_targets.py"), "--wells", w,
             "--guardband-grid", "0", "--ipc-enabled",
             "--out-calls", out[0], "--out-mosquito", out[1],
             "--out-sensitivity", out[2], "--out-ipc-summary", out[3]],
            capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        m = pd.read_csv(out[1]).set_index("sample_code")
        assert bool(m.loc["S1", "Pf_eligible"])
        assert not bool(m.loc["S2", "Pf_eligible"])
        assert m.loc["S2", "Pf_call"] == "indeterminate"


def test_ipc_enabled_without_ipc_wells_aborts():
    """Silently treating a missing IPC as a pass would be the one failure
    mode worse than not having an IPC at all."""
    wells = pd.DataFrame({
        "role": ["sample"] * 3, "plate_class": ["routine"] * 3,
        "plate": ["R01"] * 3, "well": ["A1", "A2", "A3"],
        "sample_code": ["S1"] * 3, "target": ["Pf", "Pm", "Po"],
        "cq": [25.0, None, None], "cq_conf": [0.9, 0.0, 0.0],
        "amp_status": ["Amp", "No Amp", "No Amp"],
    })
    with tempfile.TemporaryDirectory() as d:
        w = os.path.join(d, "wells.csv")
        wells.to_csv(w, index=False)
        out = [os.path.join(d, f) for f in ("c.csv", "m.csv", "s.csv")]
        r = subprocess.run(
            [sys.executable, str(BIN / "call_targets.py"), "--wells", w,
             "--guardband-grid", "0", "--ipc-enabled",
             "--out-calls", out[0], "--out-mosquito", out[1],
             "--out-sensitivity", out[2]],
            capture_output=True, text=True)
        assert r.returncode != 0
        assert "no wells carry Target Name" in r.stderr


def test_wide_table_column_names_are_explicit():
    """The wide table must be built by name, not by unpacking a pandas
    MultiIndex whose level order varies across versions. A reversed order
    silently produced 'call_Pf' instead of 'Pf_call'."""
    wells = pd.DataFrame({
        "role": ["sample"] * 6,
        "plate_class": ["routine"] * 6,
        "plate": ["R01"] * 6,
        "well": ["A1", "A2", "A3", "B1", "B2", "B3"],
        "sample_code": ["S1", "S1", "S1", "S2", "S2", "S2"],
        "target": ["Pf", "Pm", "Po"] * 2,
        "cq": [25.0, 30.0, 31.0, None, None, None],
        "cq_conf": [0.9] * 6,
        "amp_status": ["Amp"] * 3 + ["No Amp"] * 3,
    })
    with tempfile.TemporaryDirectory() as d:
        w = os.path.join(d, "wells.csv")
        wells.to_csv(w, index=False)
        out = [os.path.join(d, f) for f in ("c.csv", "m.csv", "s.csv")]
        r = subprocess.run(
            [sys.executable, str(BIN / "call_targets.py"), "--wells", w,
             "--out-calls", out[0], "--out-mosquito", out[1],
             "--out-sensitivity", out[2],
             # single-value grid: this test checks column naming, not the
             # guardband, and the no-op guard would otherwise fire on a
             # fixture this small
             "--guardband-grid", "0"],
            capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        m = pd.read_csv(out[1])
        for t in ("Pf", "Pm", "Po"):
            assert f"{t}_call" in m.columns, f"missing {t}_call: {list(m.columns)}"
            assert f"{t}_cq" in m.columns
        assert m.loc[m.sample_code == "S1", "Pf_call"].iloc[0] == "positive"
        # Po co-detected with Pf is excluded by default
        assert m.loc[m.sample_code == "S1", "Po_call"].iloc[0] == "excluded_codetection"


def test_validation_plates_are_excluded_from_prevalence():
    wells = pd.DataFrame({
        "role": ["sample"] * 6,
        "plate_class": ["routine"] * 3 + ["validation"] * 3,
        "plate": ["R01"] * 3 + ["V01"] * 3,
        "well": ["A1", "A2", "A3", "A1", "A2", "A3"],
        "sample_code": ["S1"] * 3 + ["V1"] * 3,
        "target": ["Pf", "Pm", "Po"] * 2,
        "cq": [25.0, None, None, 20.0, None, None],
        "cq_conf": [0.9] * 6,
        "amp_status": ["Amp", "No Amp", "No Amp"] * 2,
    })
    with tempfile.TemporaryDirectory() as d:
        w = os.path.join(d, "wells.csv")
        wells.to_csv(w, index=False)
        out = [os.path.join(d, f) for f in ("c.csv", "m.csv", "s.csv")]
        r = subprocess.run(
            [sys.executable, str(BIN / "call_targets.py"), "--wells", w,
             "--guardband-grid", "0", "--out-calls", out[0],
             "--out-mosquito", out[1], "--out-sensitivity", out[2]],
            capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        m = pd.read_csv(out[1])
        assert m.sample_code.tolist() == ["S1"]


def test_metadata_key_pattern_extracts_composite_ids():
    """Plate exports carry NS_003; some metadata carries NS_003_NIH00124L."""
    s = pd.Series(["NS_003_NIH00124L", "CS_001_AP005", "FT_400_XYZ"])
    assert normalise_key(s).tolist() == ["NS003NIH00124L", "CS001AP005", "FT400XYZ"]
    out = normalise_key(s, r"^([A-Za-z]+_?\d+)").tolist()
    assert out == ["NS003", "CS001", "FT400"]
    # a plain key is unaffected by the same pattern
    assert normalise_key(pd.Series(["NS_003"]), r"^([A-Za-z]+_?\d+)").tolist() == ["NS003"]


def test_ingest_drift_guard_is_wired():
    r = subprocess.run([sys.executable, str(BIN / "ingest_plate.py"), "--help"],
                       capture_output=True, text=True)
    assert "--allow-ingest-drift" in r.stdout
