# 18S qPCR TaqMan Plate Plan

## Plasmodium Detection in Blood-Fed *Anopheles gambiae* s.l.

### MIQE-Compliant Design

**Assay type:** Singleplex TaqMan qPCR (three separate reactions per sample — one per target)
**Targets:** *Plasmodium falciparum* (FAM), *P. malariae* (VIC), *P. ovale* (NED)
**Template:** Whole-mosquito gDNA diluted to 2 ng/µL; 2 µL per reaction (≈ 4 ng total)
**Plate format:** 96-well, fast or standard block

---

## 1. MIQE Compliance Summary

| MIQE Requirement                            | Status                                 | Implementation                                              |
| ------------------------------------------- | -------------------------------------- | ----------------------------------------------------------- |
| Multi-point standard curve per target       | Included                               | 6-point, 10-fold serial dilution (STD1–STD6)               |
| Standards in ≥ duplicate                   | Included                               | Duplicate (routine) / Triplicate (validation)               |
| No-template control (NTC) per target        | Included                               | 2 wells per target per plate                                |
| Extraction blank per extraction batch       | Included                               | 2 wells per target per plate                                |
| Positive control per target                 | Included                               | 2 wells per target per plate                                |
| Uninfected host negative control per target | Included                               | 2 wells per target per plate                                |
| Mosquito internal positive control (IPC)    | **Planned — not yet validated** | See §10 — development plan and plate integration strategy |
| Plate and extraction batch recording        | Required                               | See §9 record-keeping template                             |
| Sample distribution across sites / zones    | Required                               | See §8 randomisation scheme                                |

---

## 2. Primers, Probes, and Amplicon Details

| Target            | Oligo               | Sequence (5′→3′)                                              | Amplicon (bp) | Reference                 |
| ----------------- | ------------------- | ---------------------------------------------------------------- | ------------- | ------------------------- |
| *P. falciparum* | Forward             | `ATTGCTTTTGAGAGGTTTTGTTACTTT`                                  | 95            | Veron*et al.*, 2009     |
|                   | Reverse             | `GCTGTAGTATTCAAACACAATGAACTCAA`                                |               |                           |
|                   | Probe (FAM-NFQ-MGB) | `CATAACAGACGGGTAGTCAT`                                         |               |                           |
| *P. malariae*   | Forward             | `AGTTAAGGGAGTGAAGACGATCAGA`                                    | 166           | Shokoples*et al.*, 2009 |
|                   | Reverse             | `CAACCCAAAGACTTTGATTTCTCATAA`                                  |               |                           |
|                   | Probe (VIC-NFQ-MGB) | `ATGAGTGTTTCTTTTAGATAGC`                                       |               |                           |
| *P. ovale*      | Forward             | `AACCCAAAGACTTTGATTTCTCAATAA`                                  | 130           | Shokoples*et al.*, 2009 |
|                   | Reverse             | `CCGACTAGGTTTTGGATGAAAGATTTTT`                                 |               |                           |
|                   | Probe (NED-NFQ-MGB) | `CGAAAGGAATTTTCTTATT`                                          |               |                           |

> **Provenance note on the *P. malariae* and *P. ovale* oligonucleotides.**
> Two cells in this table previously held two sequences each, separated by a
> line break. The *P. malariae* reverse cell held
> `CAACCCAAAGACTTTGATTTCTCATAA` twice, which is an unambiguous duplication
> and has been collapsed. The *P. ovale* forward cell held
> `AACCCAAAGACTTTGATTTCTCAATAA` followed by a second copy of the
> *P. malariae* reverse sequence; the stray copy has been removed and the
> *P. ovale*-specific sequence retained.
>
> **What this does not settle.** The two retained sequences differ at only a
> few positions, which is what species-discriminating primers at a conserved
> site look like — but with the duplicated text removed, no single sequence
> now appears in both the *P. malariae* and the *P. ovale* row. The
> co-detection exclusion rule the analysis applies (README, *The calling
> rule*) rests on these two assays sharing a conserved reverse primer, so
> the forward/reverse assignment for both species must be verified against
> Shokoples *et al.* (2009) Table 1 and against thesis Table 3.3 before
> publication, and this table corrected to match. The sequences are
> reproduced here as the laboratory record holds them and have not been
> inferred.

---

## 3. Reaction Setup

### 3.1 Master Mix Composition (per 12 µL reaction)

A separate master mix is prepared for each target. Volumes are per well; scale by the number of wells plus 10% overage.

| Reagent                                | Stock conc. | Volume (µL)    | Final conc.   |
| -------------------------------------- | ----------- | --------------- | ------------- |
| Luna® Universal Probe qPCR Master Mix | 2×         | 6.00            | 1×           |
| Target-specific forward primer         | 10 µM      | 0.40            | 0.4 µM       |
| Target-specific reverse primer         | 10 µM      | 0.40            | 0.4 µM       |
| Target-specific probe                  | 10 µM      | 0.50            | 0.5 µM       |
| Nuclease-free water                    | —          | 2.70            | —            |
| **Sub-total (master mix)**       |             | **10.00** |               |
| Sample DNA (added separately)          | 2 ng/µL    | 2.00            | ≈ 4 ng total |
| **Total reaction volume**        |             | **12.00** |               |

> NTC wells: substitute 2 µL nuclease-free water for sample DNA.
> Extraction blank wells: use 2 µL of the blank eluate instead of sample DNA.

### 3.2 Thermocycler Programme

| Stage                                    | Step | Temp (°C) | Time     | Cycles |
| ---------------------------------------- | ---- | ---------- | -------- | ------ |
| Reverse transcriptase inactivation / UDG | 1    | 50         | 00:02:00 | 1      |
| Initial denaturation                     | 1    | 95         | 00:02:00 | 1      |
| PCR — denaturation                      | 1    | 95         | 00:00:03 | 45     |
| PCR — annealing/extension               | 2    | 60         | —       | 45     |

**Cover temperature:** 105 °C | **Reaction volume:** 12 µL

---

## 4. Standard Curve

A mixed-plasmid standard (all three targets combined) is serially diluted 10-fold in nuclease-free water containing 10 ng/µL salmon-sperm or sheared lambda DNA as a carrier to reduce adsorption losses at low concentrations.

| Standard ID | Conc. (ng/µL) | Genome Copies / Reaction (2 µL input) | Log₁₀ |
| ----------- | -------------- | -------------------------------------- | ------- |
| STD1        | 1              | 3.92 × 10⁷                           | 7.59    |
| STD2        | 0.1            | 3.92 × 10⁶                           | 6.59    |
| STD3        | 0.01           | 3.92 × 10⁵                           | 5.59    |
| STD4        | 0.001          | 3.92 × 10⁴                           | 4.59    |
| STD5        | 0.0001         | 3.92 × 10³                           | 3.59    |
| STD6        | 1 × 10⁻⁵    | 3.92 × 10²                           | 2.59    |
| (STD7)*     | 1 × 10⁻⁶    | ~39                                    | ~1.6    |
| (STD8)*     | 1 × 10⁻⁷    | ~4                                     | ~0.6    |

> \* STD7 and STD8 are optional; include on validation plates to characterise the lower limit of detection (LLOD). Omit on routine plates to maximise sample throughput. Acceptance criteria: R² ≥ 0.99, efficiency 90–110%, slope −3.58 to −3.10.

---

## 5. Controls and Their Purpose

| Control                                  | Abbreviation             | Purpose                                                                           | Min. replicates                                 |
| ---------------------------------------- | ------------------------ | --------------------------------------------------------------------------------- | ----------------------------------------------- |
| No-template control                      | NTC                      | Detects reagent contamination in the PCR setup                                    | 2 (triplicate preferred during validation)      |
| Extraction blank                         | ExB                      | Detects contamination introduced during DNA extraction                            | 2 per extraction batch represented on the plate |
| Species positive control                 | Pos_Pf / Pos_Pm / Pos_Po | Confirms assay functionality for each target; flags inhibition or reagent failure | 2                                               |
| Uninfected*Anopheles* negative control | Neg_Ctl                  | Confirms absence of cross-reactivity with the mosquito genome                     | 2                                               |
| Mosquito internal positive control (IPC) | IPC                      | Confirms successful extraction and absence of PCR inhibition at sample level      | Planned: 2 per plate once validated             |

> The same uninfected *Anopheles* gDNA stock is used as the Neg_Ctl across all three target assays on a plate. Label the stock, aliquot at −80 °C, and record the lot used per plate.

---

## 6. Plate Layouts

### Layout A — Routine Analytical Plate

**Standards in duplicate | 12 samples per plate (tested for all three targets)**

> This is the layout all 60 routine plates in this study were run on, and the
> one the sample sheet's `assignment` tab records, as a slot per specimen
> and the three wells that slot occupies. The pipeline checks every row's
> wells against the mapping below, so a plate set up from a mistyped layout
> is caught before any analysis.

Each species occupies a 4-column block (32 wells). Blocks are arranged left-to-right: *P. falciparum* (cols 1–4), *P. malariae* (cols 5–8), *P. ovale* (cols 9–12).

#### Full 96-well map

|     Row     | **1** | **2** | **3** | **4** | **5** | **6** | **7** | **8** | **9** | **10** | **11** | **12** |
| :---------: | :---------: | :---------: | :---------: | :---------: | :---------: | :---------: | :---------: | :---------: | :---------: | :----------: | :----------: | :----------: |
| **A** |   Pf_STD1   |   Pf_STD1   |   NTC_Pf   |   Spl_05   |   Pm_STD1   |   Pm_STD1   |   NTC_Pm   |   Spl_05   |   Po_STD1   |   Po_STD1   |    NTC_Po    |    Spl_05    |
| **B** |   Pf_STD2   |   Pf_STD2   |   NTC_Pf   |   Spl_06   |   Pm_STD2   |   Pm_STD2   |   NTC_Pm   |   Spl_06   |   Po_STD2   |   Po_STD2   |    NTC_Po    |    Spl_06    |
| **C** |   Pf_STD3   |   Pf_STD3   |   ExB_Pf   |   Spl_07   |   Pm_STD3   |   Pm_STD3   |   ExB_Pm   |   Spl_07   |   Po_STD3   |   Po_STD3   |    ExB_Po    |    Spl_07    |
| **D** |   Pf_STD4   |   Pf_STD4   |   ExB_Pf   |   Spl_08   |   Pm_STD4   |   Pm_STD4   |   ExB_Pm   |   Spl_08   |   Po_STD4   |   Po_STD4   |    ExB_Po    |    Spl_08    |
| **E** |   Pf_STD5   |   Pf_STD5   |   Spl_01   |   Spl_09   |   Pm_STD5   |   Pm_STD5   |   Spl_01   |   Spl_09   |   Po_STD5   |   Po_STD5   |    Spl_01    |    Spl_09    |
| **F** |   Pf_STD6   |   Pf_STD6   |   Spl_02   |   Spl_10   |   Pm_STD6   |   Pm_STD6   |   Spl_02   |   Spl_10   |   Po_STD6   |   Po_STD6   |    Spl_02    |    Spl_10    |
| **G** |   Pos_Pf   |   Neg_Ctl   |   Spl_03   |   Spl_11   |   Pos_Pm   |   Neg_Ctl   |   Spl_03   |   Spl_11   |   Pos_Po   |   Neg_Ctl   |    Spl_03    |    Spl_11    |
| **H** |   Pos_Pf   |   Neg_Ctl   |   Spl_04   |   Spl_12   |   Pos_Pm   |   Neg_Ctl   |   Spl_04   |   Spl_12   |   Pos_Po   |   Neg_Ctl   |    Spl_04    |    Spl_12    |

#### Per-block well count (Layout A)

| Content                               | Wells per species |
| ------------------------------------- | :---------------: |
| Standards (6 pts × 2 reps)           |        12        |
| NTC                                   |         2         |
| Extraction blank                      |         2         |
| Positive control                      |         2         |
| Uninfected Anopheles negative control |         2         |
| **Analytical samples**          |   **12**   |
| **Total**                       |   **32**   |

---

### Layout B — Validation / Early Optimisation Plate

**Standards in triplicate | 5 samples per plate (tested for all three targets)**

Use this layout during initial assay validation and optimisation runs. Same 4-column block structure; the third standard column is repurposed to carry the triplicate, with controls and samples distributed across rows G–H and the bottom of col 4.

#### Full 96-well map

|     Row     | **1** | **2** | **3** | **4** | **5** | **6** | **7** | **8** | **9** | **10** | **11** | **12** |
| :---------: | :---------: | :---------: | :---------: | :---------: | :---------: | :---------: | :---------: | :---------: | :---------: | :----------: | :----------: | :----------: |
| **A** |   Pf_STD1   |   Pf_STD1   |   Pf_STD1   |   NTC_Pf   |   Pm_STD1   |   Pm_STD1   |   Pm_STD1   |   NTC_Pm   |   Po_STD1   |   Po_STD1   |   Po_STD1   |    NTC_Po    |
| **B** |   Pf_STD2   |   Pf_STD2   |   Pf_STD2   |   NTC_Pf   |   Pm_STD2   |   Pm_STD2   |   Pm_STD2   |   NTC_Pm   |   Po_STD2   |   Po_STD2   |   Po_STD2   |    NTC_Po    |
| **C** |   Pf_STD3   |   Pf_STD3   |   Pf_STD3   |   NTC_Pf   |   Pm_STD3   |   Pm_STD3   |   Pm_STD3   |   NTC_Pm   |   Po_STD3   |   Po_STD3   |   Po_STD3   |    NTC_Po    |
| **D** |   Pf_STD4   |   Pf_STD4   |   Pf_STD4   |   Pos_Pf   |   Pm_STD4   |   Pm_STD4   |   Pm_STD4   |   Pos_Pm   |   Po_STD4   |   Po_STD4   |   Po_STD4   |    Pos_Po    |
| **E** |   Pf_STD5   |   Pf_STD5   |   Pf_STD5   |   Pos_Pf   |   Pm_STD5   |   Pm_STD5   |   Pm_STD5   |   Pos_Pm   |   Po_STD5   |   Po_STD5   |   Po_STD5   |    Pos_Po    |
| **F** |   Pf_STD6   |   Pf_STD6   |   Pf_STD6   |   Neg_Ctl   |   Pm_STD6   |   Pm_STD6   |   Pm_STD6   |   Neg_Ctl   |   Po_STD6   |   Po_STD6   |   Po_STD6   |   Neg_Ctl   |
| **G** |   ExB_Pf   |   ExB_Pf   |   Neg_Ctl   |   Spl_01   |   ExB_Pm   |   ExB_Pm   |   Neg_Ctl   |   Spl_01   |   ExB_Po   |    ExB_Po    |   Neg_Ctl   |    Spl_01    |
| **H** |   Spl_02   |   Spl_03   |   Spl_04   |   Spl_05   |   Spl_02   |   Spl_03   |   Spl_04   |   Spl_05   |   Spl_02   |    Spl_03    |    Spl_04    |    Spl_05    |

#### Per-block well count (Layout B)

| Content                                           | Wells per species | Total (3 species) |
| ------------------------------------------------- | :---------------: | :---------------: |
| Standards (6 pts × 3 reps)                       |        18        |        54        |
| NTC (triplicate)                                  |         3         |         9         |
| Extraction blank (duplicate)                      |         2         |         6         |
| Positive control (duplicate)                      |         2         |         6         |
| Uninfected Anopheles negative control (duplicate) |         2         |         6         |
| **Analytical samples**                      |    **5**    |   **15**   |
| **Total**                                   |   **32**   |   **96**   |

> **Note on validation samples:** Use well-characterised specimens — e.g., a dilution series of a known-positive gDNA (to assess sensitivity), a WHO/NEQAS proficiency panel, or samples with confirmed microscopy or RDT results. Spl_01–05 should represent known positives and negatives that bracket the expected working range of the assay.

> **Alternative for more validation samples:** If single-target validation is preferred (e.g., to optimise Pf independently before multiplexing), dedicate the full plate to one target and use a 6-column standard block (triplicate) + 6 columns of samples, giving 48 sample wells per plate.

---

## 7. Differences from Current (Incomplete) Plate Plan

The existing `plate_plan.xlsx` accommodates 24 samples per plate but is not fully MIQE-compliant. The table below documents what was missing and how each gap is addressed.

| Gap in existing plan                    | MIQE requirement     | Resolution in this plan                                          |
| --------------------------------------- | -------------------- | ---------------------------------------------------------------- |
| Standards run in singlet                | ≥ duplicate         | Standards in duplicate (Layout A) or triplicate (Layout B)       |
| No positive control per target          | Required             | Pos_Pf, Pos_Pm, Pos_Po added (rows G–H, col 1 per block)        |
| No uninfected mosquito negative control | Required             | Neg_Ctl added (rows G–H, col 2 per block)                       |
| No extraction blank                     | Required             | ExB added (rows C–D, col 3 per block)                           |
| No mosquito IPC                         | Strongly recommended | Noted as a current limitation; to be developed                   |
| Sample throughput (24 → 12 per plate)  | —                   | Throughput reduction is necessary to meet MIQE; plan accordingly |

---

## 8. Sample Distribution and Randomisation

### Objective

Ensure that plate number is not confounded with collection site or ecological zone so that any plate-level technical variation (efficiency, baseline, degraded reagent) is not mistaken for a biological signal.

### Preferred approach — Complete randomisation

1. Assign each sample a random number (e.g., `sample(nrow(df))` in R).
2. Sort by random number and assign samples to plates sequentially (Sample 1 → plate 1, slot 1; Sample 2 → plate 1, slot 2; … Sample 13 → plate 2, slot 1; etc.).
3. Verify post-assignment that no single plate is dominated by one site or ecological zone using a cross-tabulation (`table(plate, site)` and `table(plate, eco_zone)`).

### Fallback — Stratified block randomisation (if complete randomisation is not possible)

Use this when, for example, samples are extracted in batches that correspond to sites.

1. Divide samples into strata defined by the crossing of site × ecological zone.
2. Within each stratum, randomly assign samples to plates so that each plate receives approximately the same number of samples from each stratum.
3. If a stratum has fewer samples than plates, ensure no single plate receives all samples from that stratum.

```r
# Example R code — stratified assignment
library(dplyr)
set.seed(42)           # for reproducibility; record seed in lab notebook

samples <- samples |>
  mutate(strat = paste(site, eco_zone, sep = "_")) |>
  group_by(strat) |>
  mutate(rand_order = sample(n())) |>
  ungroup() |>
  arrange(strat, rand_order) |>
  mutate(plate_id = ((row_number() - 1) %/% 12) + 1,   # 12 samples per plate
         slot_id  =  ((row_number() - 1) %% 12) + 1)

# Verify balance
table(samples$plate_id, samples$site)
table(samples$plate_id, samples$eco_zone)
```

> Record the random seed, R version, and date of randomisation in the laboratory notebook and commit the assignment table to version control before plates are set up.

---

## 9. Record-Keeping

The following fields must be recorded for every plate to enable plate-effect QC and sensitivity analyses in the statistical model.

### Per-plate metadata (to be filled at the bench)

```
Plate ID           : ___________
Date run           : ___________
Operator           : ___________
qPCR instrument    : ___________
Block type         : Fast / Standard
Master mix lot #   : ___________
Primer / probe lots: Pf-F _______ | Pf-R _______ | Pf-probe _______
                     Pm-F _______ | Pm-R _______ | Pm-probe _______
                     Po-F _______ | Po-R _______ | Po-probe _______
Standard lot / prep date : ___________
Extraction batches on this plate : ___________
Extraction blank ID(s)           : ___________
Positive control ID / lot        : ___________
Uninfected Anopheles control lot : ___________
Notes                            : ___________
```

### Standard curve acceptance criteria (record per plate)

| Metric                      | Acceptance threshold                            | Plate value |
| --------------------------- | ----------------------------------------------- | ----------- |
| R²                         | ≥ 0.99                                         | ___         |
| Slope                       | −3.58 to −3.10 (90–110% efficiency)          | ___         |
| Efficiency                  | 90–110%                                        | ___         |
| NTC                         | No amplification, or Ct > 40                    | ___         |
| Pos control (all 3 targets) | Amplification within ± 2 Ct of historical mean | ___         |
| Neg control                 | No amplification                                | ___         |
| Extraction blank            | No amplification, or Ct > 40                    | ___         |

If any acceptance criterion is not met, the plate must be flagged, re-run where possible, and the failure documented. Do not exclude wells without documentation.

---

## 10. Notes and Limitations

### Mosquito internal positive control (IPC)

A mosquito IPC amplifies a constitutively expressed single-copy region of the *Anopheles* genome in every sample well. It simultaneously confirms (a) that gDNA extraction was successful and (b) that no PCR inhibitor co-purified with the sample. Without an IPC, a well with no Ct cannot be called a true negative — it may represent a failed extraction or inhibited reaction.

**Current status:** The IPC is not yet validated. All negative results from this assay must be reported with the caveat that failed extractions cannot be excluded. This is a named limitation for publication and must appear in the Methods section.

#### 10.1 Recommended target gene

| Candidate             | Gene ID (VectorBase) | Rationale                                                                                                           |
| --------------------- | -------------------- | ------------------------------------------------------------------------------------------------------------------- |
| *Actin-1* (primary) | AGAP006782           | Most widely used reference gene in*An. gambiae* qPCR; constitutive, single-copy, validated across sibling species |
| *RPS7*              | AGAP010592           | Alternative if*Actin-1* primers show any background in Plasmodium controls                                        |
| *RpL19*             | AGAP001397           | Second alternative; stable expression documented in blood-fed females                                               |

Design the amplicon to span an intron if possible (80-150 bp product) to discriminate genomic DNA from any residual RNA. Confirm by BLAST that the primer pair does not match *P. falciparum* 3D7, *P. malariae* UG01, or *P. ovale* curtisi/wallikeri genomes (PlasmoDB / NCBI).

#### 10.2 Fluorophore and channel

The existing assay occupies FAM (Pf), VIC (Pm), and NED (Po). ROX is reserved as the passive reference normaliser by Luna Universal Probe Master Mix and must not be used as a reporter.

| Instrument channels                                  | Recommended IPC reporter            | Notes                                                                                                                                   |
| ---------------------------------------------------- | ----------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| 5-channel (e.g. QuantStudio 5/7 Pro with Cy5 filter) | Cy5                                 | Run IPC in the same well as the Plasmodium assay. Add Cy5-labelled probe + IPC primers to each sample master mix at 0.4/0.4/0.5 uM      |
| 4-channel only (FAM/VIC/NED/ROX)                     | Run as separate 4th reaction series | Set up an additional plate section or separate plate using a FAM-labelled IPC probe in a singleplex reaction without Plasmodium primers |

#### 10.3 IPC primer/probe design parameters

| Parameter             | Target value                                       |
| --------------------- | -------------------------------------------------- |
| Amplicon size         | 80-150 bp                                          |
| Primer length         | 18-22 nt                                           |
| Primer Tm             | 58-62 deg C                                        |
| Probe Tm              | >= 68 deg C (>= 10 deg above primer Tm)            |
| Probe label           | Cy5-NFQ-MGB (or FAM-NFQ-MGB for separate reaction) |
| GC content            | 40-60%                                             |
| No 3-prime G on probe | Required (reduces background ligation artefacts)   |

#### 10.4 IPC plate integration (once validated)

When using a Cy5 IPC in the same well as the Plasmodium assay:

- **Standard wells (STD1-6):** No IPC required. Plasmid standards do not contain Anopheles genomic sequence; adding IPC would not be interpretable.
- **NTC wells:** IPC must be **negative** (no template = no IPC amplification). If IPC amplifies in NTC, reagent contamination is indicated.
- **Extraction blank wells:** IPC must be **negative**. If IPC amplifies, the extraction reagents or consumables are contaminated with Anopheles gDNA.
- **Neg_Ctl (uninfected Anopheles) wells:** IPC must **amplify**. Failure indicates poor-quality gDNA in the negative control stock.
- **Pos_Ctl wells:** IPC must **amplify**. Failure indicates inhibition or poor extraction of the positive control stock.
- **Sample wells:** IPC must **amplify** (Ct <= defined threshold). If IPC is absent or Ct > threshold, the sample is flagged as extraction failure or inhibited and must not be reported as Plasmodium-negative. Resample or re-extract if material is available.

The Layout A plate map (Section 6) will be updated to include IPC designation in each sample and control well once the IPC assay is validated and the fluorophore channel confirmed.

#### 10.5 IPC development and validation plan

| Step | Action                                                                           | Acceptance criterion                                        |
| ---- | -------------------------------------------------------------------------------- | ----------------------------------------------------------- |
| 1    | Design Actin-1 primers and Cy5 probe; order oligos                               | Primer-BLAST: no off-target in Plasmodium genomes           |
| 2    | Singleplex validation of IPC alone (serial dilution of An. gambiae gDNA)         | R2 >= 0.99, efficiency 90-110%                              |
| 3    | Multiplex compatibility: run IPC + each of Pf, Pm, Po assays pairwise            | Ct shift <= 0.5 cycles vs. singleplex; no probe cross-talk  |
| 4    | Establish pass/fail Ct threshold on 20+ field extractions                        | Threshold set at mean Ct + 2 SD from extraction-only blanks |
| 5    | Run 10 extraction blanks with IPC                                                | All blanks: IPC negative (undetermined)                     |
| 6    | Validate interpretation rules using archived samples with known infection status | Sensitivity and specificity unchanged vs. assay without IPC |

#### 10.6 Reporting until IPC is validated

Include the following statement in the Methods section of any manuscript:

> "An internal positive control targeting the *Anopheles gambiae Actin-1* locus (AGAP006782) is under development but was not incorporated during the study period. Consequently, samples yielding no amplification for all three Plasmodium targets cannot be definitively classified as true negatives, as failed extractions cannot be excluded. The number of samples returning undetermined Ct values for all targets is reported in Table S[X] to allow readers to assess the potential impact of this limitation."

### Singleplex design and throughput

Running three separate singleplex reactions per sample (rather than a triplex) improves analytical sensitivity and avoids probe-to-probe competition but reduces throughput from 24 samples (original plan) to 12 samples per plate. For the full sample set, plan accordingly:

```
Number of plates needed = ceil(N_samples / 12)
```

Plan for 10–15% re-run capacity on top of this to accommodate plate failures and ambiguous results.

### Extraction batch confounding

If it is operationally impossible to achieve complete randomisation (e.g., samples are extracted by site-specific field teams), plate number and extraction batch must be entered as random effects or covariates in the statistical model. Coordinate with the bioinformatics / statistical analysis plan to ensure these fields are consistently captured in the sample metadata file.

---

## Abbreviations

| Abbreviation | Meaning                                                                                                |
| ------------ | ------------------------------------------------------------------------------------------------------ |
| ExB          | Extraction blank                                                                                       |
| IPC          | Internal positive control                                                                              |
| LLOD         | Lower limit of detection                                                                               |
| MGB          | Minor groove binder (quencher moiety)                                                                  |
| MIQE         | Minimum Information for Publication of Quantitative Real-Time PCR Experiments (Bustin*et al.*, 2009) |
| Neg_Ctl      | Uninfected*Anopheles gambiae* s.l. negative control                                                  |
| NFQ          | Non-fluorescent quencher                                                                               |
| NTC          | No-template control                                                                                    |
| Pf / Pm / Po | *Plasmodium falciparum* / *P. malariae* / *P. ovale*                                             |
| Pos_Pf/Pm/Po | Species-specific positive control                                                                      |
| Spl          | Analytical sample (unknown)                                                                            |
| STD          | Standard (plasmid dilution series)                                                                     |

---

*References:*
Bustin SA *et al.* (2009) The MIQE guidelines: minimum information for publication of quantitative real-time PCR experiments. *Clin Chem* 55(4):611–622.
Shokoples SE *et al.* (2009) Multiplexed real-time PCR for discriminating *Plasmodium* species with improved sensitivity over other real-time and conventional nested PCR-based methods. *J Clin Microbiol* 47(4):975–980.
Veron V *et al.* (2009) Real-time PCR detection of *Plasmodium* in human and *Anopheles* samples. *Malaria J* (as cited in plate plan metadata).
