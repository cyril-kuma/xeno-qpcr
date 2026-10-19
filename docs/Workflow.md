# Workflow reference

One page per stage: what it reads, what it writes, and what it will refuse to
let through. The README covers the reasoning; this is the mechanical
reference.

Exit codes are distinct so a failure is identifiable from a log line alone.

| code | stage | meaning |
|---|---|---|
| 2 | `INGEST_PLATE`, `MERGE_WELLS` | ingest drift: raw and ingested amplification-status counts disagree |
| 3 | `MERGE_WELLS` | two exports resolve to the same plate identifier |
| 3 | `CALL_TARGETS` | the guardband is a no-op for a non-dominant target |
| 4 | `IMPORT_SAMPLE_SHEET`, `CALL_TARGETS`, `RANDOMISE` | an input is unusable: a missing column, a well that contradicts its slot, two specimens in one well, more unmatched specimens than `--max-unmatched`, a specimen join with zero overlap, an allocation exceeding plate capacity |
| 5 | `INGEST_PLATE` | role integrity: `Sample Name` and `Task` disagree, or a well resolves to no role |
| 6 | `VALIDATE_PLATE_MAPS` | the plates that ran do not match the assignment |

---

## RANDOMISE — `bin/randomise.py`

Design-time. Reached only with `--randomise`.

**Reads** the specimen metadata, and optionally a first-round plate layout
(`--first_round_slots`) with `Sample Name` and `plate` columns.

**Writes** `sample_assignment_new.csv` and `allocation_balance.csv`.

**Algorithm.** Strata are zone × site × year × collection method. Specimens
are shuffled within each stratum; strata are ordered largest first with ties
broken at random; specimens are dealt across plates in laps with a fresh
random plate permutation each lap. Spare slots become inter-plate replicates,
drawn one each from the largest strata and placed on a plate other than the
specimen's own. Slot order within each plate is then randomised. Each slot
maps to three wells by the Layout A block structure.

**Refuses** an allocation whose specimen count exceeds `n_plates ×
slots_per_plate`, naming the number of plates needed.

**Warns** when a stratifying variable is absent from the metadata, because the
allocation is then unverified on it.

**Note.** This does not reproduce the frozen allocation bit-for-bit unless the
inputs and the NumPy generator match those used at the time. The allocation
this study ran on lives in the sample sheet's `assignment` tab.

## IMPORT_SAMPLE_SHEET — `bin/import_sample_sheet.py`

**Reads** `data/sample_sheet.xlsx`: a `specimens` tab and an `assignment`
tab. Finds each header row by looking for `sample_code` or `plate_id` in the
first four rows, because a one-line note sits above it.

**Writes** `specimens.csv`, `sample_assignment.csv` (one row per slot, with
`is_replicate`, `replicate_of_plate` and `specimen_status` derived),
`plate_map_expected.csv` (long: one row per plate × well × target, with the
specimen the design expects there), `specimen_reconciliation.csv`.

**Checks**
- *slot against wells*: each row's well triplet against the Layout A mapping
  its `slot_id` implies. The two are typed independently, so one contradicts
  the other and the error is a transcription slip in the sheet.
- *specimens against assignment*, in both directions: specimens run without a
  record, and specimens described but never run. Both counted and named.

**Refuses**
- either tab missing a required column, listing what was found
- `sample_code` not unique in the specimens tab — it is the join key everywhere
- non-numeric `slot_id` or `run_order`
- a well that contradicts its slot
- an allocation placing two specimens in the same well
- more assigned specimens without a record than `--max-unmatched` (default
  20). A handful is a documented gap; a large number means the two tabs
  describe different specimens, and a silent intersection could not tell the
  two apart

**Does not drop** specimens with no record. They are physically in wells on
real plates, so the plate-map gate must see them; removing them here would
make their wells look unexpected. They leave the analysis at `CALL_TARGETS`
under `restrict_to_metadata`.

**Does not read** the laboratory's own tracking workbook. The run log, QC
checklist, reagent log and manual results entry are bench records with no
analytical role, and hand-transcribed results would be a second source of
truth for numbers the instrument exports already carry.

## INGEST_PLATE — `bin/ingest_plate.py`

One process per plate.

**Reads** one QuantStudio `.xlsx`. Tries calamine first, then openpyxl: 57 of
the 60 exports here carry a malformed `workbookView` attribute that openpyxl
rejects outright. Finds the header row by scanning for `Well` and
`Target Name`, falling back to `--header-row`.

**Writes** `wells_<plate>.csv` and `audit_<plate>.json`.

**Derives** `plate` from the file name, keeping a repeat suffix
(`Plas_18S_qPCR_R07-2.xlsx` → `R07-2`); `plate_class` from whether the plate
identifier matches `--validation-plate-pattern`; `role` from `Sample Name`
(see `bin/well_roles.py`); `control_target` from a control's suffix;
`sample_code` from the specimen name, normalised and left null on controls.

**Refuses**
- a file name with no plate identifier at the end
- role integrity failure: a well whose `Sample Name` matches no control
  pattern while `Task` marks it a control, or whose implied role is not one
  its `Task` can carry. The offending wells are printed.
- ingest drift, unless `--allow-ingest-drift`

## MERGE_WELLS — `bin/merge_ingest.py`

**Writes** `wells_merged.csv` and `ingest_audit.csv` (one row per plate).

**Refuses** duplicate plate identifiers, and any plate with ingest drift.

## VALIDATE_PLATE_MAPS — `bin/validate_plate_maps.py`

The gate. Its output is the well table every later stage reads, so the
dataflow enforces the order.

**Reads** `wells_merged.csv` and `plate_map_expected.csv`.

**Writes** `wells.csv` and `plate_map_validation.csv`.

**Checks** every specimen well on every routine plate against the assignment
on (plate, well, target), classifying disagreements as `misplaced`,
`wrong_slot_same_plate`, `unexpected` or `absent`. Plates in the assignment
with no export are reported as `plates_designed_not_run` rather than as
mismatches. Validation plates are exempt.

**Refuses** any `misplaced`, `unexpected` or `absent` well, unless
`--allow-plate-map-mismatch`.

## PLATE_QC — `bin/plate_qc.py`

**Writes** `target_qc.csv` (one row per plate × target, with a separate
failure count per control class and a `fail_reasons` string) and
`control_qc.csv` (one row per control well).

**Gates** each control class separately — see the README table. Positive
controls additionally warn when more than `--pos-ctl-cq-tolerance` cycles from
the run-wide per-target median; a warning is not a failure.

**Reports** control wells whose dedicated target disagrees with the well's
`Target Name`, which would mean a control was pipetted into the wrong block.

**Acceptance** R² ≥ `--r2-min` (0.99), efficiency within `--eff-min`/`--eff-max`,
slope within `--slope-min`/`--slope-max`. `--ctl-amplification-max` decides
whether a control amplified and is independent of the calling rule's
`--cq-trace-max`.

## CROSS_REACTIVITY — `bin/cross_reactivity.py`

**Reads** `wells.csv`, using only wells on validation plates.

**Writes** `cross_reactivity.json` (the verdict) and `cross_reactivity.csv`
(the per-check detail).

**Verdict** `clean`, `cross_reactive` or `not_assessed`. Evidence considered:
uninfected-host gDNA amplifying in any *Plasmodium* channel, and a positive
control for one species amplifying in another species' channel.

`not_assessed` is not a failure. It records that specificity rests on the
per-plate uninfected-host negative controls and the published primer design
rather than on a validation plate in this run.

## CALL_TARGETS — `bin/call_targets.py`

**Reads** `wells.csv`, `specimens.csv`, and the cross-reactivity verdict.

**Writes** `target_calls.csv` (well level), `mosquito_level.csv` (specimen
level, joined to the metadata), `crossamp_sensitivity.csv`, `ipc_summary.csv`.

**Order of operations** matters and is fixed: base calls → the pre-exclusion
snapshot → IPC routing (if enabled) → co-detection exclusion → eligibility
and detection → the guardband arm, computed on the pre-exclusion snapshot.

**Excludes** wells on validation plates before any call is made.

**Refuses**
- a well table missing a required column
- no wells with `role == 'sample'`
- a requested target with no wells at all
- a specimen join with zero overlap, after scanning every column of the
  specimens tab with and without pattern extraction and naming the flags that
  would work
- `--ipc-enabled` with no wells carrying the IPC target
- a guardband that is a no-op for any non-dominant target that has
  co-detections to act on

## REPLICATE_CONCORDANCE — `bin/replicate_concordance.py`

**Reads** `target_calls.csv` and `sample_assignment.csv`; optionally a
first-round mosquito-level table (`--round1-calls`).

**Writes** `replicate_pairs.csv` and `replicate_summary.csv`.

**Independence** is checked before agreement is counted: a pair whose Cq *and*
Cq-confidence agree to within `--independence-eps` is a duplicated record,
not a second measurement, and is excluded with status `not_independent`. A
pair that returned no amplification on both plates is uninformative rather
than duplicated. `call_agreement` and the Cq differences are reported as
`not estimable` when no comparison has two independent Cq values.

## RANDOMISATION_CHECK — `bin/randomisation_check.py`

**Reads** `mosquito_level.csv`, `wells.csv`, `sample_assignment.csv` and
`specimens.csv`. The assignment records position and identity only, so the
attributes the allocation was balanced on are joined from the specimens on
`sample_code`. The balance is tested over every slot as allocated, including
specimens with no record: restricting to the analysed set first would test
the balance of a subset rather than of the randomisation.

**Writes** `allocation_balance.csv` (one row per stratifying variable),
`batch_diagnostic.csv`, `plate_positivity.csv`.

**Tests** plate against each variable in `--strata` by chi-square, with the
largest share any level takes on any plate; then counts plates at exactly 0%
or 100% positivity against a permutation null over plates of the same sizes;
then reports the one-way ANOVA intraclass correlation of detection for each
variable in `--icc-groups`.

**Flag** `CONFOUNDED` (the allocation itself is imbalanced), `PLATE_EFFECT`
(balanced allocation, clustered detection — a technical effect, so carry plate
as a random effect or cluster variable), or `ok`. The `interpretation` row
states which, in words, and `make_report.py` reproduces it.

## REPORT — `bin/make_report.py`

**Writes** `prevalence_by_zone.csv` and `analysis_report.md`.

Prevalence is per zone and overall, for each target, for any *Plasmodium* and
for mixed detection, with Wilson intervals. The report carries run integrity,
the per-class control counts, the cross-reactivity verdict, prevalence, the
sensitivity arm, the allocation balance, the residual plate structure with its
interpretation, the replicate concordance, and the IPC counts.
