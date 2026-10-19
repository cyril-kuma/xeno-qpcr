#!/usr/bin/env nextflow
/*
 * xeno-qpcr — Plasmodium 18S singleplex qPCR detection in blood-fed
 * Anopheles gambiae sensu lato.
 *
 * The workflow follows the order the study was carried out in:
 *
 *   randomisation  ->  sample sheet  ->  plate maps  ->  plate exports
 *                  ->  QC  ->  calls  ->  checks  ->  report
 *
 * Two inputs describe the whole experiment: the plate exports, and one
 * sample sheet holding the specimens and the slot each was run in. The
 * laboratory's own tracking workbook stays with the bench.
 *
 * The first stage is design-time and runs once, before any plate is set up.
 * It is reached only with `--randomise`, so an analysis run can never
 * silently re-allocate specimens; an analysis starts instead from the
 * allocation the plates were actually set up from.
 *
 * The strict syntax has no `-entry`, so the mode is a parameter and the entry
 * workflow dispatches on it.
 *
 * Stage order is enforced by the dataflow rather than by convention:
 * VALIDATE_PLATE_MAPS emits the well table every later stage reads, so
 * nothing downstream can run on plates that do not match the assignment.
 */
nextflow.enable.dsl = 2

include { RANDOMISE; IMPORT_SAMPLE_SHEET        } from './modules/design.nf'
include { INGEST_PLATE; MERGE_WELLS;
          VALIDATE_PLATE_MAPS                   } from './modules/ingest.nf'
include { PLATE_QC; CROSS_REACTIVITY            } from './modules/qc.nf'
include { CALL_TARGETS; REPLICATE_CONCORDANCE;
          RANDOMISATION_CHECK; REPORT           } from './modules/analysis.nf'

def helpMessage() {
    return """
    xeno-qpcr ${workflow.manifest.version}

    Analysis of an executed design
      nextflow run main.nf -params-file conf/params.yaml -profile conda

      --plates        glob for QuantStudio .xlsx exports (quote it)
      --sample_sheet  workbook with a 'specimens' tab and an 'assignment' tab
      --round1_calls  optional first-round mosquito-level table, for the
                      wider replicate concordance comparison
      --outdir        ${params.outdir}

    Design-time allocation (run once, before any plate is set up)
      nextflow run main.nf --randomise \\
          --specimens data/sample_sheet.xlsx \\
          --randomisation_seed 20261001

      --specimens           a CSV, or a sample sheet whose 'specimens' tab
                            holds the mosquitoes to allocate
      --n_plates            ${params.n_plates}
      --slots_per_plate     ${params.slots_per_plate}
      --first_round_slots   optional first-round layout, supplying
                            first_run_plate

    Calling rule
      --cq_positive_max     ${params.cq_positive_max}   Cq at or below which a target is positive
      --cq_trace_max        ${params.cq_trace_max}      upper bound of the trace band
      --cqconf_targets      ${params.cqconf_targets}    targets that must also clear the Cq-confidence floor
      --cqconf_min          ${params.cqconf_min}        that floor
      --codetection_exclude ${params.codetection_exclude}      cannot be positive when ${params.dominant_target} is
      --guardband_grid      ${params.guardband_grid}    cross-amplification sensitivity arm
      --ipc_enabled         ${params.ipc_enabled}    route IPC failures to indeterminate

    Standard-curve and control acceptance
      --miqe_r2_min            ${params.miqe_r2_min}
      --miqe_eff_min/max       ${params.miqe_eff_min} / ${params.miqe_eff_max}
      --miqe_slope_min/max     ${params.miqe_slope_min} / ${params.miqe_slope_max}
      --ctl_amplification_max  ${params.ctl_amplification_max}     a control with a Cq at or below this amplified
      --pos_ctl_cq_tolerance   ${params.pos_ctl_cq_tolerance}      cycles either side of the run median

    Sample-sheet reconciliation
      --max_unmatched_specimens  ${params.max_unmatched_specimens}    assigned specimens with no
                                 record before the run is stopped
    """.stripIndent()
}

/*
 * Design-time allocation. Reached only with `--randomise`: its output is a
 * proposal to be reviewed, committed and written into a sample sheet, not an
 * input the analysis picks up on its own. RANDOMISE prints the notebook and
 * version-control reminders on its own stdout.
 */
workflow ALLOCATE {
    main:
    if( !params.specimens ) { error('--specimens is required with --randomise') }

    slots_ch = params.first_round_slots
               ? Channel.fromPath(params.first_round_slots, checkIfExists: true)
               : Channel.fromPath("${projectDir}/assets/NO_FIRST_ROUND")

    RANDOMISE(file(params.specimens), slots_ch)
}

/*
 * Analysis of an executed design.
 */
workflow ANALYSE {
    main:
    if( !params.plates )       { error('--plates is required') }
    if( !params.sample_sheet ) { error('--sample_sheet is required: the workbook holding the specimens and the slot each was run in') }

    // ---- 1. the specimens, and where each was run ----
    IMPORT_SAMPLE_SHEET(file(params.sample_sheet))

    // ---- 2. the plate exports ----
    plates_ch = Channel.fromPath(params.plates, checkIfExists: true)
                       .map { f -> tuple(f.baseName.tokenize('_').last(), f) }

    INGEST_PLATE(plates_ch)
    MERGE_WELLS(INGEST_PLATE.out.wells.collect(),
                INGEST_PLATE.out.audit.collect())

    // ---- 3. the gate: the plates that ran are the plates designed ----
    //
    // Validated against the FULL assignment, including specimens with no
    // record in the specimens tab. Those specimens are physically in wells on
    // real plates, so dropping them before this point would make their wells
    // look unexpected and fail the gate. They leave the analysis at
    // CALL_TARGETS instead, under restrict_to_metadata.
    VALIDATE_PLATE_MAPS(MERGE_WELLS.out.wells, IMPORT_SAMPLE_SHEET.out.plate_map)
    wells_ch = VALIDATE_PLATE_MAPS.out.wells

    // ---- 4. QC and the specificity evidence ----
    PLATE_QC(wells_ch)
    CROSS_REACTIVITY(wells_ch)

    // ---- 5. calls ----
    CALL_TARGETS(wells_ch, IMPORT_SAMPLE_SHEET.out.specimens,
                 CROSS_REACTIVITY.out.verdict)

    // ---- 6. the checks that qualify them ----
    round1_ch = params.round1_calls
                ? Channel.fromPath(params.round1_calls, checkIfExists: true)
                : Channel.fromPath("${projectDir}/assets/NO_ROUND1")

    REPLICATE_CONCORDANCE(CALL_TARGETS.out.calls,
                          IMPORT_SAMPLE_SHEET.out.assignment, round1_ch)
    RANDOMISATION_CHECK(CALL_TARGETS.out.mosquito, wells_ch,
                        IMPORT_SAMPLE_SHEET.out.assignment,
                        IMPORT_SAMPLE_SHEET.out.specimens)

    // ---- 7. report ----
    REPORT(CALL_TARGETS.out.mosquito, CALL_TARGETS.out.sensitivity,
           PLATE_QC.out.target_qc, PLATE_QC.out.control_qc,
           RANDOMISATION_CHECK.out.diagnostic,
           RANDOMISATION_CHECK.out.balance,
           MERGE_WELLS.out.audit, VALIDATE_PLATE_MAPS.out.report,
           REPLICATE_CONCORDANCE.out.summary, CALL_TARGETS.out.ipc,
           CROSS_REACTIVITY.out.detail,
           IMPORT_SAMPLE_SHEET.out.reconciliation)
}

/*
 * The entry workflow. The strict syntax has no `-entry`, so it dispatches on
 * `--randomise` rather than Nextflow selecting a named workflow.
 */
workflow {
    main:
    if( params.help ) {
        log.info(helpMessage())
        return
    }

    if( params.randomise ) {
        ALLOCATE()
    }
    else {
        ANALYSE()
    }

    // Both outcomes are reported from onComplete. An `onError:` section is
    // invoked with an argument while the section body compiles to a zero-arg
    // closure, which aborts the handler itself with "Invalid method
    // invocation `doCall` with arguments: null" and masks the task error that
    // triggered it.
    onComplete:
    log.info(workflow.success
        ? (params.randomise
            ? "Allocation written to ${params.outdir}/randomise. Review it, record the seed in the laboratory notebook, and write it into a sample sheet before setting up any plate."
            : "Completed. Results in ${params.outdir}")
        : "Failed. See .nextflow.log and the failing task's work directory.")
}
