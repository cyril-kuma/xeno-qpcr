# Changelog

## 1.1.0

The pipeline now reads one `data/sample_sheet.xlsx` — a `specimens` tab and an
`assignment` tab — in place of `data/metadata/Metadata.csv` and the
11-sheet laboratory tracking workbook. The tracking workbook stays with the
bench: it records what was done, it keeps being edited after an analysis is
frozen, and its hand-transcribed `Results_Entry` sheet was a second source of
truth for numbers the instrument exports already carry.

**Breaking.** `--tracker` and `--metadata` are replaced by `--sample_sheet`,
and `--randomise` now takes `--specimens` (a CSV, or a sample sheet whose
`specimens` tab holds them). Results are unchanged: every number in the 1.0.0
run is reproduced exactly, including the allocation chi-square statistics.

**Two checks replace one that would have become circular.** The 1.0.0 import
cross-checked the workbook's `Plate_Maps` sheet against `Sample_Assignment`,
two independent transcriptions of the same allocation. Generated from a single
source those tabs would validate nothing while still reporting a pass, so the
plate-map tab is gone and `slot_id` is instead checked against the three well
columns, which are typed independently and related by the Layout A mapping.
The specimens and assignment tabs are reconciled in both directions and the
counts reported, with `--max_unmatched_specimens` stopping a run where the
join has broken; `specimen_status` is derived from that reconciliation rather
than hand-maintained. `Balance_Check` is dropped — it was only printed, and
the balance it showed is recomputed from the assignment.

## 1.0.0

First release. The workflow takes raw QuantStudio exports of singleplex 18S
rRNA TaqMan qPCR assays for *Plasmodium falciparum*, *P. malariae* and
*P. ovale* through MIQE quality control to zone-stratified prevalence
estimates, in the order the study was carried out: randomisation, sample
assignment, plate maps, then the analysis of what was run.

**Stages.** `RANDOMISE` (design-time, reached with `--randomise`),
`IMPORT_SAMPLE_SHEET`, `INGEST_PLATE`, `MERGE_WELLS`, `VALIDATE_PLATE_MAPS`,
`PLATE_QC`, `CROSS_REACTIVITY`, `CALL_TARGETS`, `REPLICATE_CONCORDANCE`,
`RANDOMISATION_CHECK`, `REPORT`. Each is a standalone CLI in `bin/`.

**The dataset this release was validated on.** 60 routine plates, 5760 wells,
714 specimens assigned and run, 701 with metadata and analysed. All 2160
specimen wells match the randomised assignment. All 180 plate × target
standard curves pass at R² ≥ 0.99. No control failure in any of the four
classes. Plate ICC of detection 0.000 against 0.755 by collection site and
0.538 by bioclimatic zone, with 0 of 60 plates at 0% or 100% positivity.

**Two findings this release surfaces rather than absorbs.**

- Nine of the eighteen inter-plate replicate comparisons carry identical Cq
  *and* identical Cq-confidence across the two plates, so the second record
  duplicates the first rather than measuring it again. The replicate arm is
  reported as not estimable and does not support a reproducibility claim for
  this dataset.
- The *P. malariae* and *P. ovale* oligonucleotide cells in
  `docs/Assay_design.md` §2 held duplicated text. The duplication is removed
  and the provenance recorded, but the forward/reverse assignment for both
  species still needs verification against Shokoples *et al.* (2009) Table 1,
  because the co-detection exclusion rule cites that table.

**Requirements.** Nextflow 25.10 or later with the strict syntax parser
(`NXF_SYNTAX_PARSER=v2` on 25.04–25.10; the default from 26.04). Python 3.11
with pandas, numpy, scipy, statsmodels, openpyxl, python-calamine and
tabulate; `environment.yml` pins them for the `conda` profile.

**Testing.** 62 pytest regression tests, and a two-plate synthetic fixture the
`test` profile runs the whole workflow on, so CI exercises every stage without
real specimen data. The fixture's generators are in `tests/fixtures/` and CI
regenerates it rather than trusting the committed copy.
