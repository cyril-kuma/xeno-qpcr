# xeno-qpcr

Detection of *Plasmodium* DNA in blood-fed *Anopheles gambiae* sensu lato by
singleplex 18S rRNA TaqMan qPCR, from raw QuantStudio exports through MIQE
quality control to prevalence estimates.

Three targets — *P. falciparum* (FAM), *P. malariae* (VIC), *P. ovale* (NED) —
run as **separate singleplex reactions**: one target per well, three wells per
specimen, with per-target standards and controls on every plate.

Two inputs describe the whole experiment: the plate exports, and one sample
sheet holding the specimens and the slot each was run in.

**Version 1.1.0.** Requires **Nextflow 25.10 or later** and the strict syntax
parser; see *Nextflow version notes*.

---

## Quick start

```bash
# the reportable run: every parameter pinned in one file
NXF_SYNTAX_PARSER=v2 nextflow run main.nf -params-file conf/params.yaml -profile conda

# the whole workflow on a two-plate synthetic fixture, no real data needed
NXF_SYNTAX_PARSER=v2 nextflow run main.nf -profile test

# the regression suite
python3 -m pytest tests -q
```

Profiles: `standard` (local), `conda`, `docker`, `slurm`, `test`.

---

## The workflow

The stages follow the order the study was carried out in. The allocation comes
first, because the whole design rests on it; the plate-map check comes before
any analysis, because it is the only place a mis-set plate is visible.

```
                    ┌─ RANDOMISE ──────── design-time, --randomise only
                    │     allocate specimens to plates and slots
                    ▼
sample_sheet.xlsx ─ IMPORT_SAMPLE_SHEET ─ the specimens, and the slot each
                    │                     was run in; reconciled in both
                    │                     directions
  data/plates/ ─── INGEST_PLATE ───────── one process per plate; re-reads
                    │                     the raw sheet and asserts the
                    │                     ingested counts match it
                   MERGE_WELLS
                    │
                   VALIDATE_PLATE_MAPS ── THE GATE: every specimen well
                    │                     must hold the specimen the
                    │                     assignment put there
         ┌──────────┼──────────┐
         ▼          ▼          ▼
     PLATE_QC  CROSS_REACTIVITY │          MIQE curves, four control
         │          │          │          classes gated separately, and
         │          └──────────┤          the specificity verdict
         │                     ▼
         │                CALL_TARGETS ─── well calls, mosquito-level
         │                     │           collapse, sensitivity arm
         │          ┌──────────┴──────────┐
         │          ▼                     ▼
         │  REPLICATE_CONCORDANCE  RANDOMISATION_CHECK
         │          │                     │
         └──────────┴─────────┬───────────┘
                              ▼
                           REPORT ──────── prevalence with Wilson
                                           intervals, plus every check
```

Stage order is enforced by the dataflow, not by convention:
`VALIDATE_PLATE_MAPS` emits the well table every later stage reads, so nothing
downstream can run on plates that do not match the assignment.

Every stage is a standalone CLI in `bin/`, so any step can be run or debugged
outside Nextflow.

---

## 1. The allocation

Plate-level technical variation is real — a reagent lot, a run, a thermal
block — and if plate is correlated with collection site or bioclimatic zone,
that variation is indistinguishable from a biological difference between
zones. Assay design §8 therefore requires specimens to be distributed so that
plate is independent of site, zone, collection method and year.

`bin/randomise.py` does the allocation: stratified systematic assignment,
dealing each stratum across plates in laps with a fresh random plate
permutation per lap, then randomising slot order within each plate so position
carries no information either. Slots left spare are filled with a second
aliquot of a specimen already placed on another plate — the inter-plate
technical replicates.

The allocation used in this study was produced at **seed 20261001** and is
frozen in the sample sheet's `assignment` tab. An analysis run never
re-allocates: the
allocator is reached only with `--randomise`, and writes to a separate file.

```bash
NXF_SYNTAX_PARSER=v2 nextflow run main.nf --randomise \
    --specimens data/sample_sheet.xlsx \
    --randomisation_seed 20261001
```

Review the output, record the seed and the date in the laboratory notebook,
commit it, and write it into a sample sheet **before** setting up any plate.

## 2. The sample sheet, and the gate

`data/sample_sheet.xlsx` is the pipeline's only description of the experiment.
Two tabs:

| tab            | grain                        | holds                                                                |
| -------------- | ---------------------------- | -------------------------------------------------------------------- |
| `specimens`  | one row per specimen (701)   | everything known about the mosquito, independent of where it was run |
| `assignment` | one row per plate slot (720) | plate, slot, the three target wells, the stratum, the role           |

They are separate on purpose. 720 slots cover 714 specimens plus 6 inter-plate
replicates, against 701 described specimens; merging the tabs would repeat
specimen attributes on the replicate rows and lose the distinction between
specimens that were *run* and specimens that are *described* — which is the
distinction every denominator in this analysis rests on.

The laboratory's own tracking workbook — run log, QC checklist, reagent log,
manual results entry — stays with the bench and is not read here. It records
what was done, it keeps being edited after an analysis is frozen, and its
hand-transcribed results would be a second source of truth for numbers the
instrument exports already carry.

`IMPORT_SAMPLE_SHEET` makes two checks.

**Slot against wells.** `slot_id` and the three well columns are typed
independently but related by the Layout A block structure, so each row's well
triplet is checked against the mapping its slot implies. A mistyped well cell
is caught here rather than surfacing later as a plate that disagrees with its
own design.

**Specimens against assignment.** The two tabs are reconciled in both
directions and the counts reported, never silently intersected:

```
720 slots -> 714 specimens run -> 13 without a record -> 701 analysed
```

A silent inner join cannot tell a handful of known gaps from a mis-keyed
column that matches half the sheet: both yield a clean run on whatever
survived. `--max_unmatched_specimens` (default 20) stops a run where the join
has plainly broken, and the report carries the chain above so the thesis can
state its denominators. In this dataset the 13 are `FT028`–`FT031` and
`FT064`–`FT072` — two consecutive runs, so whole collection batches are
absent from the specimens tab rather than scattered rows.

Those 13 are **reported, not dropped here.** They are physically in wells on
real plates, so the plate-map gate must see them; removing them first would
make their wells look unexpected. They leave the analysis at `CALL_TARGETS`,
under `restrict_to_metadata`.

`specimen_status` is derived from this reconciliation rather than
hand-maintained in the sheet, so it cannot go stale.

`VALIDATE_PLATE_MAPS` then matches every specimen well in the exports against
the assignment on (plate, well, target). **This is the only stage that can
detect a plate set up from the wrong row.** Calls from a mis-set plate are
internally consistent and attached to the wrong specimen, so the zone, site
and collection method that specimen contributes to would all be wrong, and no
QC gate, no curve and no permutation test would notice. Four classes of
disagreement abort:

| class          | meaning                                                                                                       |
| -------------- | ------------------------------------------------------------------------------------------------------------- |
| `misplaced`  | the well holds a different specimen than the assignment expects                                               |
| `wrong_slot` | the specimen is on its assigned plate but not in its slot — a transcription slip rather than a swapped plate |
| `unexpected` | the export has a specimen well the assignment does not list                                                   |
| `absent`     | the assignment lists a well the export does not contain                                                       |

On the 60 plates shipped here, all 2160 specimen wells match.

---

## 3. Roles and controls

The instrument's `Task` column cannot separate the control classes the design
distinguishes. On every plate in this study `Task == "NTC"` covers both the
no-template controls and the extraction blanks:

| `Task`             | `Sample Name`              | role assigned        | wells/plate |
| -------------------- | ---------------------------- | -------------------- | ----------- |
| `STANDARD`         | `{Pf,Pm,Po}_STD_1`…`_6` | `standard`         | 36          |
| `UNKNOWN`          | specimen codes               | `sample`           | 36          |
| `NTC`              | `NTC_{Pf,Pm,Po}`           | `ntc`              | 6           |
| `NTC`              | `ExB_{Pf,Pm,Po}`           | `extraction_blank` | 6           |
| `Positive control` | `Pos_{Pf,Pm,Po}`           | `positive`         | 6           |
| `NEGATIVE`         | `Neg_Ctl`                  | `negative_host`    | 6           |

Roles therefore come from `Sample Name`, which is unambiguous, and `Task` is
used only to cross-check the assignment; a disagreement aborts ingest. A
consequence worth having: role assignment no longer depends on the spelling of
`Task`, so an export writing `POSITIVE CONTROL` or `Positive Control`
validates identically.

Each class is gated separately, because each implicates a different failure
and a different corrective action:

| role                 | expectation                                     | verdict on breach                  |
| -------------------- | ----------------------------------------------- | ---------------------------------- |
| `ntc`              | no amplification                                | `FAIL_pcr_reagent_contamination` |
| `extraction_blank` | no amplification                                | `FAIL_extraction_contamination`  |
| `negative_host`    | no amplification                                | `FAIL_host_cross_reactivity`     |
| `positive`         | amplification                                   | `FAIL_no_amplification`          |
| `positive`         | within ±2 Cq of the run-wide per-target median | `WARN_pos_ctl_cq_drift`          |

A single pooled "negative control failure" count reports a contaminated
extraction, contaminated PCR reagents and a primer set cross-reacting with the
mosquito genome as the same number.

**Standard-curve acceptance** follows assay design §4 and §9: R² ≥ **0.99**,
efficiency 90–110%, slope −3.58 to −3.10. All 180 plate × target curves in
this dataset pass at 0.99.

`ctl_amplification_max` (whether a control amplified) is a separate parameter
from `cq_trace_max` (the upper bound of the trace band for specimen wells).
Both are 40, but they answer different questions and moving one must not move
the other.

---

## 4. The calling rule

```
amp_status == "Inconclusive"              → invalid
amp_status == "No Amp"  or Cq is null     → negative
Cq ≤ cq_positive_max (38)                 → positive
    …except targets in cqconf_targets (Pm), which must also have
      Cq Conf ≥ cqconf_min (0.50), else → trace
cq_positive_max < Cq ≤ cq_trace_max (40)  → trace
Cq > cq_trace_max                         → negative
```

**Trace is not positive.** It is kept as its own state so the evaluable
denominator stays honest rather than silently absorbing marginal wells.

**Amplification status is authoritative.** Some exports carry a Cq on a well
the instrument called `No Amp` — plate R39 well A4 reads Cq 22.71 with
`Amp Status = No Amp`. A Cq alone is not a detection, and the rule above
treats the status as decisive. This is deliberate and locked by a test.

**Eligibility**

| denominator      | definition                                                   |
| ---------------- | ------------------------------------------------------------ |
| per-target       | `call` is neither `invalid` nor `indeterminate`        |
| `eligible_any` | at least one target positive, **or** all targets valid |
| `eligible_all` | all targets valid — the mixed-detection denominator         |

**Co-detection exclusion — the primary specificity rule.** Targets listed in
`codetection_exclude` (default `Po`) cannot be called positive in a mosquito
where the dominant target (Pf) is positive. Every such co-detection is treated
as cross-amplification rather than as mixed infection.

The justification is structural and empirical. The assay follows the Shokoples
design — species-specific forward primers with a **shared conserved reverse
primer** for Pm and Po — so discrimination rests on the forward primer and a
short MGB probe. Empirically 92.6% of Po positives co-occur with Pf, at a
median delay of +3.1 cycles: the signature of the Po reaction amplifying Pf
template. Note this is *stricter* than a guardband of 0, which retains
co-detections whose Cq happens to fall below the dominant target's — exclusion
gives 29 Po positives, a guardband of 0 gives 69.

The oligonucleotide table this rule cites needs one correction before
publication; see the provenance note in `docs/Assay_design.md` §2.

**Cross-reactivity verdict.** The rule's structural premise is checkable on a
validation plate, where known-composition material makes a signal in the wrong
channel unexplainable any other way. `CROSS_REACTIVITY` produces that verdict
rather than assuming it — `clean`, `cross_reactive`, or `not_assessed` — and
`CALL_TARGETS` carries it into the output so the report states the evidence
state instead of leaving the reader to assume a validation was done. No
validation plate is shipped here, so the verdict is `not_assessed`: specificity
rests on the per-plate uninfected-host negative controls and on the published
primer design. Drop `V01.xlsx`, `V02.xlsx` into `data/plates/` and the stage
assesses them (`validation_plate_pattern`, default `^V`); such plates are
excluded from prevalence.

**Guardband sensitivity arm.** Asks what would happen if a Cq-delay guardband
were used *instead of* exclusion, and is therefore computed on the
**pre-exclusion** calls. Running it after exclusion leaves it nothing to act on
and yields a flat column — which is what happened in an earlier analysis. The
arm never alters the primary call, and a guard aborts if any non-dominant
target's count is constant across the grid.

---

## 5. The internal positive control

Without an IPC, a well with no Cq cannot be distinguished from a failed
extraction or an inhibited reaction, so it cannot be called a true negative.
The IPC was not incorporated during the study period (assay design §10), so
`ipc_enabled` is **false** and the report instead carries the count assay
design §10.6 requires — specimens returning no signal at every target — so a
reader can bound the limitation. In this dataset that is **165 of 701
specimens (23.5%)**.

The machinery is in place for when an IPC exists. With `--ipc_enabled`, a
specimen whose IPC is absent or later than `ipc_cq_max` has its non-positive
calls become `indeterminate` and drops out of every denominator; it is never
reported as negative. Positive calls stand — a target that amplified did so
regardless of what the IPC shows. Enabling the IPC when no well carries the
IPC target aborts, because silently treating a missing IPC as a pass would be
the one failure mode worse than not having an IPC at all.

---

## 6. Did the randomisation work?

`RANDOMISATION_CHECK` answers two questions, and keeping them apart is the
point of the stage.

**Decorrelation** is a property of the allocation alone, checkable before a
plate is run: a chi-square test of independence between plate and each
stratifying variable, plus the largest share any one level takes on any one
plate. On this allocation all four variables are balanced.

**Residual plate structure** is then a separate question. Given a verified
allocation, detection that still clusters by plate cannot be zone confounding —
the allocation rules that out — so it is a technical plate effect, and the
conclusion is that plate belongs in the model as a random effect or cluster
variable rather than that prevalence is uninterpretable. Reporting the
intraclass correlation of detection by plate next to the ICC by site and by
zone makes that concrete: a plate ICC far above the zone ICC is technical, the
reverse is biology.

On this dataset:

| level            | ICC of detection |
| ---------------- | ---------------- |
| plate            | **0.000**  |
| collection site  | 0.755            |
| bioclimatic zone | 0.538            |

with 0 of 60 plates at 0% or 100% positivity against a permutation median of
1, so `batch_effect_flag` reads `ok`. The variance sits at the site and zone
level, which is where biology would put it.

`batch_effect_flag` takes three values: `CONFOUNDED` (the allocation itself is
imbalanced — resolve it before interpreting any plate-level result),
`PLATE_EFFECT` (balanced allocation, clustered detection — carry plate as a
random effect), and `ok`.

---

## 7. Inter-plate replicates

Six specimens were aliquoted onto a second plate. These pairs are the only
measurement in the design that separates plate-to-plate technical variation
from specimen-to-specimen biological variation: the same DNA, the same assay,
two plates.

Each pair is checked for **independence** before its agreement is counted. Two
reactions in two wells on two plates do not produce the same Cq to six decimal
places, and they do not produce the same Cq-confidence either; when they do,
the second record is a copy rather than a second measurement, and the pair
says nothing about reproducibility however well it agrees.

**In this dataset 9 of 18 target comparisons fail that check** — the two
plates carry identical Cq *and* identical Cq-confidence (for example FT254 Pf
reads 30.891409 on both R39 E4 and R02 D4, with Cq Conf 0.960875 on both).
The remaining 9 returned no amplification on both plates, which agrees by
construction. So `call_agreement` and `mean_abs_delta_cq` are reported as
**not estimable**, and the replicate arm does not support a reproducibility
claim for this dataset. Agreement computed over copies would have read as
perfect concordance and would have been the most misleading number in the
report.

Specimens also carry `first_run_plate`, the plate they sat on in the first
round of this assay. Supplying that round's mosquito-level table with
`--round1_calls` adds a second, larger concordance comparison across the whole
set; without it the comparison is reported as not performed, since the
first-round result table is not part of this repository.

---

## 8. Joining plates to the specimens tab

Plate exports carry the bare Sample Name (`NS_003`); the pipeline normalises it
by stripping underscores and upper-casing (`NS003`). The column named by
`--metadata_key` is normalised the same way.

Here `sample_code` joins directly and is the default. The `sample_id` column is
a composite (`NS003NIH00124L`) and would need a regex whose group 1 is the
joinable part:

```bash
--metadata_key sample_id --metadata_key_pattern '^([A-Za-z]+_?\d+)'
```

`CALL_TARGETS` reports the link count before doing any work, and on a zero
overlap it scans every metadata column, with and without pattern extraction,
and names the column and pattern that would actually join. Read the error; it
gives you the flags.

Of the 714 specimens on the plates, 701 have a record in the specimens tab and
are analysed; the other 13 are excluded by `restrict_to_metadata`, having
already been named by the reconciliation. Every per-target denominator is 701
minus that target's invalid wells, so the counts close arithmetically against
the assignment.

---

## 9. Guards you should not remove

Each exists because the fault occurred and nothing noticed.

**Ingest integrity.** An earlier implementation silently lost the NED/Po
channel on 42 of 60 plates: wells the instrument had called `Amp` arrived as
`No Amp`, with a null Cq and a Cq-confidence of exactly 0.000. Reported
*P. ovale* prevalence fell from ~60% to 4.3%, and nothing flagged it.
`ingest_plate.py` re-reads the raw sheet and asserts, per target, that
ingested amplification-status counts equal the raw counts. Drift aborts.
`--allow_ingest_drift` downgrades this to a warning **for debugging only**.

**Guardband no-op.** In an earlier analysis Po returned an identical count at
guardbands 0, 1, 2 and 3. The cause was ordering: co-detection exclusion had
already been applied upstream, so the guardband had nothing left to act on. A
threshold rule that gives the same answer at every threshold is not running.

**Denominator integrity.** An earlier implementation applied co-detection
exclusion by overwriting `amp_status` at ingest. That also converted 16
`Inconclusive` Po wells into `No Amp`, moving them from *invalid* to *negative*
and inflating the Po denominator. Applying the rule at the call level, as here,
leaves invalid wells invalid.

**Plate identity.** `merge_ingest.py` aborts if two exports resolve to the same
plate identifier, which would pool their wells. The identifier regex keeps a
repeat suffix, so `Plas_18S_qPCR_R07-2.xlsx` is plate `R07-2`, not plate `2`.

Two overrides exist for diagnosis and must never be set for a reportable run:
`--allow_ingest_drift` and `--allow_plate_map_mismatch`. Both are `false` in
`conf/params.yaml`.

---

## 10. Interpreting the output

Because extracts are whole-mosquito, a positive means *Plasmodium* DNA was
present in the mosquito. It does not distinguish parasites ingested with the
blood meal from an established infection, and **it is not a sporozoite rate**.

Nor is it a parasitaemia: the standards are plasmid dilutions, so a Cq maps to
plasmid copies in the reaction, not to parasites in the host blood meal of
unknown and unmeasured volume.

---

## Outputs

```
results/
├── import_sample_sheet/    specimens.csv, sample_assignment.csv,
│                           plate_map_expected.csv,
│                           specimen_reconciliation.csv
├── validate_plate_maps/    wells.csv, plate_map_validation.csv
├── plate_qc/               target_qc.csv, control_qc.csv
├── cross_reactivity/       cross_reactivity.json, cross_reactivity.csv
├── call_targets/           target_calls.csv, mosquito_level.csv,
│                           crossamp_sensitivity.csv, ipc_summary.csv
├── replicate_concordance/  replicate_pairs.csv, replicate_summary.csv
├── randomisation_check/    batch_diagnostic.csv, allocation_balance.csv,
│                           plate_positivity.csv
├── report/                 analysis_report.md, prevalence_by_zone.csv
└── pipeline_info/          report.html, timeline.html, trace.txt, dag.html
```

`results/randomise/` holds `sample_assignment_new.csv` and
`allocation_balance.csv` when run with `--randomise`.

---

## Repository layout

```
xeno-qpcr/
├── main.nf                 workflow; dispatches on --randomise
├── nextflow.config         params, profiles, resources
├── conf/params.yaml        pinned parameters for the reportable run
├── modules/                design.nf · ingest.nf · qc.nf · analysis.nf
├── bin/                    standalone CLIs, also usable without Nextflow
│   ├── well_roles.py            the role vocabulary, shared by ingest and QC
│   ├── randomise.py             design-time allocation
│   ├── import_sample_sheet.py   sample sheet → flat CSVs, reconciled
│   ├── ingest_plate.py          one QuantStudio export → tidy wells
│   ├── merge_ingest.py          concatenate + flatten audits
│   ├── validate_plate_maps.py   the gate
│   ├── plate_qc.py              MIQE curves, four control classes
│   ├── cross_reactivity.py      specificity verdict from validation plates
│   ├── call_targets.py          calls, collapse, sensitivity arm, IPC
│   ├── replicate_concordance.py inter-plate pairs, with independence check
│   ├── randomisation_check.py   decorrelation and residual plate structure
│   └── make_report.py           prevalence with Wilson intervals
├── data/sample_sheet.xlsx  specimens + assignment, the only description
│                           of the experiment the pipeline reads
├── data/plates/            60 QuantStudio exports (R01–R60)
├── docs/Assay_design.md    the MIQE plate plan this pipeline implements
├── tests/                  pytest suite (62 tests)
│   ├── fixtures/                generators for the synthetic test data
│   └── data/                    the two-plate fixture the test profile runs
├── assets/                 placeholders standing in for optional inputs
├── environment.yml
├── CITATION.cff · CHANGELOG.md · LICENSE
└── .github/workflows/ci.yml
```

---

## Before publishing

- Verify the *P. malariae* / *P. ovale* oligonucleotide orientations against
  Shokoples *et al.* (2009) Table 1 and correct `docs/Assay_design.md` §2.
  The co-detection rule cites that table.
- Resolve the nine non-independent replicate comparisons against the
  instrument run files, or state the limitation in the Methods.
- Decide whether `data/plates/` and `data/sample_sheet.xlsx` ship with the
  repository.
  They are included here because the analysis is only reproducible with them;
  `.gitignore` has commented-out lines to exclude them if the specimen data
  must stay private.
- Replace the GitHub owner in `nextflow.config` (`docker` profile container),
  `CITATION.cff` and the manifest `homePage` if the repository moves.
- Tag `v1.1.0` and archive via Zenodo for a DOI.

---

## Nextflow version notes

This workflow targets the **strict syntax**, optional from 25.04 and the
default from 26.04. On 25.04–25.10, set `NXF_SYNTAX_PARSER=v2`; from 26.04 it
is the default and the variable is unnecessary.

Four constructs are worth knowing if you adapt the code.

**Event handlers.** `workflow.onComplete { ... }` at the top level of a script
is a *statement*, and the strict syntax forbids mixing statements with
declarations. Since 25.10 the handler is an `onComplete:` section inside the
entry workflow:

```groovy
workflow {
    main:
    // ...

    onComplete:
    log.info("done")
}
```

Using `main:` becomes mandatory once the workflow has any other section, and
**only the entry workflow may carry one** — a named workflow with an
`onComplete:` section fails to compile.

**No `onError:` section.** Nextflow invokes the handler with one argument
while an `onError:` section body compiles to a zero-arg closure, so the
handler aborts with `Invalid method invocation doCall with arguments: null` —
and because it runs during failure handling, it masks the task error that
triggered it. Report both outcomes from `onComplete:` using `workflow.success`,
as this workflow does.

**No `-entry`.** The strict syntax removes it; a named workflow is selected by
a parameter from the entry workflow instead, which is why the allocator is
reached with `--randomise` rather than `-entry randomise`.

**`exit` is gone.** Use `error('message')` to abort with a message, and a bare
`return` to stop the entry workflow cleanly, as `--help` does here.

One configuration behaviour to know: a `--param` given on the command line
**replaces a profile's whole `params` block** rather than merging into it. The
`test` profile therefore sets its own `outdir`, and should be run on its own or
overridden with `-params-file`, not with `--param` flags.
