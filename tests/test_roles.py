"""Role assignment, and the control classes the Task column cannot separate."""
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

from well_roles import (  # noqa: E402
    AMPLIFICATION_EXPECTED, CONTROL_ROLES, EXTRACTION_BLANK, NEGATIVE_HOST,
    NO_AMPLIFICATION_EXPECTED, NTC, POSITIVE, SAMPLE, STANDARD, assign_role,
    control_target, role_from_sample_name)


@pytest.mark.parametrize("name,task,role", [
    # On every plate in this study Task == "NTC" covers both the no-template
    # controls and the extraction blanks. Sample Name is what separates them.
    ("NTC_Pf", "NTC", NTC),
    ("NTC_Pm", "NTC", NTC),
    ("NTC_Po", "NTC", NTC),
    ("ExB_Pf", "NTC", EXTRACTION_BLANK),
    ("ExB_Pm", "NTC", EXTRACTION_BLANK),
    ("ExB_Po", "NTC", EXTRACTION_BLANK),
    ("Pos_Pf", "Positive control", POSITIVE),
    ("Neg_Ctl", "NEGATIVE", NEGATIVE_HOST),
    ("Pf_STD_1", "STANDARD", STANDARD),
    ("NS_147", "UNKNOWN", SAMPLE),
])
def test_roles_come_from_sample_name(name, task, role):
    assert assign_role(name, task) == (role, "ok")


def test_ntc_and_extraction_blank_are_distinct_roles():
    """The whole point of assigning from Sample Name. An NTC that amplifies
    implicates the PCR reagents; an extraction blank that amplifies
    implicates the extraction. Collapsing them loses that distinction."""
    ntc, _ = assign_role("NTC_Pf", "NTC")
    exb, _ = assign_role("ExB_Pf", "NTC")
    assert ntc != exb
    assert {ntc, exb} <= set(CONTROL_ROLES)
    assert ntc in NO_AMPLIFICATION_EXPECTED
    assert exb in NO_AMPLIFICATION_EXPECTED
    assert NO_AMPLIFICATION_EXPECTED[ntc] != NO_AMPLIFICATION_EXPECTED[exb]


def test_every_no_amplification_class_has_its_own_verdict():
    verdicts = list(NO_AMPLIFICATION_EXPECTED.values())
    assert len(verdicts) == len(set(verdicts)) == 3
    assert set(NO_AMPLIFICATION_EXPECTED) == {NTC, EXTRACTION_BLANK,
                                              NEGATIVE_HOST}
    assert set(AMPLIFICATION_EXPECTED) == {POSITIVE}


@pytest.mark.parametrize("task", [
    "Positive control", "POSITIVE CONTROL", "positive control",
    "Positive Control", "  Positive control  ",
])
def test_task_casing_does_not_affect_the_role(task):
    """The old mapping keyed on an upper-cased Task against the literal
    "POSITIVE CONTROL", so an export writing "Positive control" fell through
    to a catch-all and the positive controls stopped being gated. Roles now
    come from Sample Name and the Task cross-check is case-folded."""
    assert assign_role("Pos_Pf", task) == (POSITIVE, "ok")


@pytest.mark.parametrize("name", ["Pf_STD1", "Pf_STD_1", "Po_STD_6",
                                  "Pm_STD6", "pf_std_3"])
def test_standard_naming_variants(name):
    """Assay design §6 writes Pf_STD1; the exports write Pf_STD_1."""
    assert role_from_sample_name(name) == STANDARD


def test_disagreement_between_sample_name_and_task_is_reported():
    role, why = assign_role("NTC_Pf", "UNKNOWN")
    assert role == NTC
    assert why != "ok" and "UNKNOWN" in why


def test_unrecognised_control_name_does_not_pass_as_a_sample():
    """A control the patterns do not know must not be silently analysed as a
    specimen, which is what a Task-only mapping would do."""
    role, why = assign_role("Blank_7", "NTC")
    assert role == "unclassified"
    assert why != "ok"


def test_control_target_comes_from_the_suffix():
    assert control_target("NTC_Pm", NTC) == "Pm"
    assert control_target("ExB_Po", EXTRACTION_BLANK) == "Po"
    assert control_target("Pos_Pf", POSITIVE) == "Pf"
    # One uninfected-host stock serves all three target assays (§5), so the
    # name carries no target and the well's Target Name is the only source.
    assert control_target("Neg_Ctl", NEGATIVE_HOST) is None
