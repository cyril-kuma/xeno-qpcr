"""Control gating, standard-curve acceptance, and replicate independence."""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))


def _wells(rows):
    base = {"plate": "R01", "plate_class": "routine", "quantity": None,
            "r2": 0.995, "slope": -3.3, "efficiency": 100.0,
            "y_intercept": 40.0, "control_target": None, "cq_conf": 0.9,
            "amp_status": "Amp", "sample_code": None, "well": "A1"}
    return pd.DataFrame([{**base, **r} for r in rows])


def _standards(target="Pf"):
    return [{"role": "standard", "target": target, "sample_name": f"{target}_STD_{i}",
             "cq": 12.0 + 3 * i, "quantity": 10.0 ** (7 - i), "well": f"A{i}"}
            for i in range(1, 7)]


def _run_qc(wells_df, extra=()):
    with tempfile.TemporaryDirectory() as d:
        w = os.path.join(d, "wells.csv")
        t = os.path.join(d, "target_qc.csv")
        c = os.path.join(d, "control_qc.csv")
        wells_df.to_csv(w, index=False)
        r = subprocess.run(
            [sys.executable, str(BIN / "plate_qc.py"), "--wells", w,
             "--out-target", t, "--out-control", c, *extra],
            capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        return r, pd.read_csv(t), pd.read_csv(c)


def test_each_control_class_gets_its_own_verdict():
    """A contaminated extraction, contaminated PCR reagents and a primer set
    cross-reacting with the host genome are three different findings. One
    pooled 'negative control failure' count reports them as the same."""
    wells = _wells(_standards() + [
        {"role": "ntc", "target": "Pf", "sample_name": "NTC_Pf",
         "control_target": "Pf", "cq": 33.0},
        {"role": "extraction_blank", "target": "Pf", "sample_name": "ExB_Pf",
         "control_target": "Pf", "cq": 35.0},
        {"role": "negative_host", "target": "Pf", "sample_name": "Neg_Ctl",
         "cq": 36.0},
        {"role": "positive", "target": "Pf", "sample_name": "Pos_Pf",
         "control_target": "Pf", "cq": None, "amp_status": "No Amp"},
    ])
    _r, tq, cq = _run_qc(wells)
    verdicts = dict(zip(cq.role, cq.verdict))
    assert verdicts["ntc"] == "FAIL_pcr_reagent_contamination"
    assert verdicts["extraction_blank"] == "FAIL_extraction_contamination"
    assert verdicts["negative_host"] == "FAIL_host_cross_reactivity"
    assert verdicts["positive"] == "FAIL_no_amplification"

    row = tq.iloc[0]
    for role in ("ntc", "extraction_blank", "negative_host", "positive"):
        assert row[f"n_fail_{role}"] == 1, f"{role} not counted separately"
    assert row.n_control_failures == 4
    assert not row.curve_ok
    for role in ("ntc", "extraction_blank", "negative_host", "positive"):
        assert role in row.fail_reasons


def test_clean_controls_pass():
    wells = _wells(_standards() + [
        {"role": "ntc", "target": "Pf", "sample_name": "NTC_Pf",
         "control_target": "Pf", "cq": None, "amp_status": "No Amp"},
        {"role": "extraction_blank", "target": "Pf", "sample_name": "ExB_Pf",
         "control_target": "Pf", "cq": None, "amp_status": "No Amp"},
        {"role": "negative_host", "target": "Pf", "sample_name": "Neg_Ctl",
         "cq": None, "amp_status": "No Amp"},
        {"role": "positive", "target": "Pf", "sample_name": "Pos_Pf",
         "control_target": "Pf", "cq": 20.0},
    ])
    _r, tq, cq = _run_qc(wells)
    assert set(cq.verdict) == {"pass"}
    assert bool(tq.iloc[0].curve_ok)


def test_a_control_above_the_amplification_threshold_has_not_amplified():
    """Assay design §9 accepts an NTC with no amplification or Ct > 40."""
    wells = _wells(_standards() + [
        {"role": "ntc", "target": "Pf", "sample_name": "NTC_Pf",
         "control_target": "Pf", "cq": 41.0},
    ])
    _r, _tq, cq = _run_qc(wells)
    assert cq.loc[cq.role == "ntc", "verdict"].iloc[0] == "pass"


def test_control_threshold_is_independent_of_the_trace_band():
    """ctl_amplification_max and cq_trace_max share a value but answer
    different questions, so moving one must not move the other."""
    wells = _wells(_standards() + [
        {"role": "ntc", "target": "Pf", "sample_name": "NTC_Pf",
         "control_target": "Pf", "cq": 39.0},
    ])
    _r, _tq, strict = _run_qc(wells, extra=("--ctl-amplification-max", "38"))
    assert strict.loc[strict.role == "ntc", "verdict"].iloc[0] == "pass"
    _r, _tq, loose = _run_qc(wells, extra=("--ctl-amplification-max", "40"))
    assert loose.loc[loose.role == "ntc", "verdict"].iloc[0].startswith("FAIL")


def test_r2_threshold_is_the_design_value():
    """Assay design §4 and §9 require R2 >= 0.99. A curve at 0.985 passed
    under the earlier 0.98 setting."""
    wells = _wells([{**s, "r2": 0.985} for s in _standards()])
    _r, tq, _cq = _run_qc(wells)
    assert not bool(tq.iloc[0].curve_ok)
    assert "r2=" in tq.iloc[0].fail_reasons
    _r, relaxed, _cq = _run_qc(wells, extra=("--r2-min", "0.98"))
    assert bool(relaxed.iloc[0].curve_ok)


def test_a_drifting_positive_control_warns_without_failing():
    wells = _wells(_standards() + [
        {"role": "positive", "target": "Pf", "sample_name": "Pos_Pf",
         "control_target": "Pf", "cq": 20.0, "well": "G1"},
        {"role": "positive", "target": "Pf", "sample_name": "Pos_Pf",
         "control_target": "Pf", "cq": 20.2, "well": "H1"},
        {"role": "positive", "target": "Pf", "sample_name": "Pos_Pf",
         "control_target": "Pf", "cq": 26.0, "well": "G2"},
    ])
    _r, tq, cq = _run_qc(wells)
    v = cq.loc[cq.role == "positive", "verdict"].tolist()
    assert "WARN_pos_ctl_cq_drift" in v
    # a warning is not a failure: the curve still passes
    assert int(tq.iloc[0].n_fail_positive) == 0
    assert bool(tq.iloc[0].curve_ok)


# --------------------------------------------------------------- replicates --

def _run_conc(calls, assignment, extra=()):
    with tempfile.TemporaryDirectory() as d:
        c = os.path.join(d, "calls.csv")
        a = os.path.join(d, "assignment.csv")
        p = os.path.join(d, "pairs.csv")
        s = os.path.join(d, "summary.csv")
        calls.to_csv(c, index=False)
        assignment.to_csv(a, index=False)
        r = subprocess.run(
            [sys.executable, str(BIN / "replicate_concordance.py"),
             "--calls", c, "--assignment", a, "--targets", "Pf",
             "--out-pairs", p, "--out-summary", s, *extra],
            capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        return r, pd.read_csv(p), pd.read_csv(s)


def _pair(primary_cq, replicate_cq, primary_conf=0.9, replicate_conf=0.9):
    calls = pd.DataFrame([
        {"plate": "R01", "well": "E3", "target": "Pf", "role": "sample",
         "sample_code": "S1", "call": "positive", "cq": primary_cq,
         "cq_conf": primary_conf},
        {"plate": "R02", "well": "F3", "target": "Pf", "role": "sample",
         "sample_code": "S1", "call": "positive", "cq": replicate_cq,
         "cq_conf": replicate_conf},
    ])
    asg = pd.DataFrame([{
        "plate_id": "R02", "slot_id": 2, "sample_code": "S1",
        "well_Pf": "F3", "well_Pm": "F7", "well_Po": "F11",
        "is_replicate": True, "replicate_of_plate": "R01",
        "stratum": "x", "first_run_plate": "P01"}])
    return calls, asg


def test_identical_cq_and_confidence_is_flagged_as_not_independent():
    """Two reactions in two wells on two plates do not give the same Cq to
    six decimals. When they do, the second record is a copy, and agreement
    computed over copies would read as perfect reproducibility."""
    _r, pairs, summary = _run_conc(*_pair(30.891409, 30.891409))
    assert pairs.status.tolist() == ["not_independent"]
    got = dict(zip(summary.statistic, summary.value))
    assert int(got["n_comparisons_not_independent"]) == 1
    assert got["call_agreement"] == "not estimable"
    assert "must not be cited" in got["replicate_arm_verdict"]


def test_two_genuine_measurements_are_compared():
    _r, pairs, summary = _run_conc(
        *_pair(30.89, 31.14, primary_conf=0.96, replicate_conf=0.94))
    assert pairs.status.tolist() == ["compared"]
    assert bool(pairs.independent.iloc[0])
    got = dict(zip(summary.statistic, summary.value))
    assert float(got["call_agreement"]) == 1.0
    assert float(got["mean_abs_delta_cq"]) == 0.25


def test_no_amplification_on_both_plates_is_uninformative_not_a_copy():
    """A pair that returned nothing on both plates agrees by construction.
    It is not evidence of a duplicated record, and it is not evidence of
    reproducibility either."""
    _r, pairs, summary = _run_conc(
        *_pair(None, None, primary_conf=0.0, replicate_conf=0.0))
    assert pairs.status.tolist() == ["compared"]
    got = dict(zip(summary.statistic, summary.value))
    assert int(got["n_comparisons_with_two_cq_values"]) == 0
    assert got["call_agreement"] == "not estimable"
    assert "by construction" in got["replicate_arm_verdict"]
