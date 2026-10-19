#!/usr/bin/env python3
"""Well role vocabulary, shared by ingest and QC.

The instrument's `Task` column carries five values, and they do not separate
the control classes the assay design distinguishes. On every plate in this
study `Task == "NTC"` covers both the no-template controls and the extraction
blanks:

    Task                Sample Name            wells/plate
    STANDARD            {Pf,Pm,Po}_STD_1..6    36
    UNKNOWN             specimen codes         36
    NTC                 NTC_{Pf,Pm,Po}          6
    NTC                 ExB_{Pf,Pm,Po}          6
    Positive control    Pos_{Pf,Pm,Po}          6
    NEGATIVE            Neg_Ctl                 6

Those two NTC classes answer different questions. An NTC that amplifies
implicates the PCR reagents; an extraction blank that amplifies implicates the
extraction reagents or consumables, and the corrective action differs. Assay
design §5 and §9 list them as separate acceptance criteria, so the pipeline
assigns roles from `Sample Name`, which is unambiguous, and then uses `Task`
only to cross-check the assignment.

Role assignment is therefore independent of the `Task` spelling. The
cross-check compares case-folded values, so an export writing
"POSITIVE CONTROL" or "Positive Control" validates identically.
"""
from __future__ import annotations

import re

SAMPLE = "sample"
STANDARD = "standard"
NTC = "ntc"
EXTRACTION_BLANK = "extraction_blank"
POSITIVE = "positive"
NEGATIVE_HOST = "negative_host"

#: Control classes that must show no amplification, each with the failure it
#: implicates. Keys are roles; values are the verdict recorded on a breach.
NO_AMPLIFICATION_EXPECTED = {
    NTC: "FAIL_pcr_reagent_contamination",
    EXTRACTION_BLANK: "FAIL_extraction_contamination",
    NEGATIVE_HOST: "FAIL_host_cross_reactivity",
}

#: Control classes that must amplify.
AMPLIFICATION_EXPECTED = {POSITIVE: "FAIL_no_amplification"}

CONTROL_ROLES = tuple(NO_AMPLIFICATION_EXPECTED) + tuple(AMPLIFICATION_EXPECTED)

# Ordered; the first match wins. Standard IDs are accepted with or without the
# separator before the level number, so both Pf_STD1 (assay design §6) and
# Pf_STD_1 (as the exports write it) resolve to `standard`.
_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"^NTC(?:[_-].*)?$", re.I), NTC),
    (re.compile(r"^ExB(?:[_-].*)?$", re.I), EXTRACTION_BLANK),
    (re.compile(r"^Pos(?:[_-].*)?$", re.I), POSITIVE),
    (re.compile(r"^Neg[_-]?Ctl(?:[_-].*)?$", re.I), NEGATIVE_HOST),
    (re.compile(r"^(?:Pf|Pm|Po)[_-]STD[_-]?\d+$", re.I), STANDARD),
]

#: Role each `Task` value is expected to carry, case-folded.
_TASK_EXPECTED = {
    "unknown": {SAMPLE},
    "standard": {STANDARD},
    "ntc": {NTC, EXTRACTION_BLANK},
    "positive control": {POSITIVE},
    "negative": {NEGATIVE_HOST},
}


def role_from_sample_name(sample_name: str) -> str | None:
    """Role implied by the Sample Name, or None when no pattern matches."""
    s = str(sample_name).strip()
    for pattern, role in _PATTERNS:
        if pattern.match(s):
            return role
    return None


def assign_role(sample_name: str, task: str) -> tuple[str, str]:
    """Return (role, consistency) for one well.

    `consistency` is "ok", or a description of how Sample Name and Task
    disagree. Callers decide whether a disagreement is fatal; ingest treats it
    as fatal, because a mislabelled control silently weakens a QC gate.
    """
    named = role_from_sample_name(sample_name)
    key = str(task).strip().casefold()
    expected = _TASK_EXPECTED.get(key)

    if named is None:
        if key == "unknown":
            return SAMPLE, "ok"
        if expected is None:
            return "unclassified", (
                f"Task {task!r} is not one of {sorted(_TASK_EXPECTED)} and "
                f"Sample Name {sample_name!r} matches no control pattern"
            )
        return "unclassified", (
            f"Sample Name {sample_name!r} matches no control pattern, but "
            f"Task {task!r} marks the well as a control"
        )

    if expected is None:
        return named, (
            f"Sample Name {sample_name!r} implies role {named!r}, but Task "
            f"{task!r} is not one of {sorted(_TASK_EXPECTED)}"
        )
    if named not in expected:
        return named, (
            f"Sample Name {sample_name!r} implies role {named!r}, which is "
            f"not among the roles Task {task!r} can carry ({sorted(expected)})"
        )
    return named, "ok"


def control_target(sample_name: str, role: str) -> str | None:
    """Target a control well is dedicated to, from its suffix.

    NTC_Pf, ExB_Pm and Pos_Po name their target; Neg_Ctl is a single
    uninfected-host gDNA stock shared across all three target assays on a
    plate (assay design §5), so it names none and the well's Target Name is
    the only source.
    """
    if role == NEGATIVE_HOST:
        return None
    m = re.search(r"[_-](Pf|Pm|Po)$", str(sample_name).strip(), re.I)
    return m.group(1).title() if m else None
