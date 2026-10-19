/*
 * Calls, the checks that qualify them, and the report.
 */

process CALL_TARGETS {
    input:
    path wells
    path specimens
    path cross_reactivity

    output:
    path 'target_calls.csv',         emit: calls
    path 'mosquito_level.csv',       emit: mosquito
    path 'crossamp_sensitivity.csv', emit: sensitivity
    path 'ipc_summary.csv',          emit: ipc

    script:
    def restrict = params.restrict_to_metadata ? '--restrict-to-metadata' : ''
    def keypat = params.metadata_key_pattern
                 ? "--metadata-key-pattern '${params.metadata_key_pattern}'" : ''
    def ipc = params.ipc_enabled ? '--ipc-enabled' : ''
    """
    call_targets.py --wells ${wells} --metadata ${specimens} \\
        --cross-reactivity ${cross_reactivity} \\
        --metadata-key '${params.metadata_key}' ${keypat} \\
        --targets '${params.targets}' --dominant '${params.dominant_target}' \\
        --cq-positive-max ${params.cq_positive_max} \\
        --cq-trace-max ${params.cq_trace_max} \\
        --cqconf-targets '${params.cqconf_targets}' \\
        --cqconf-min ${params.cqconf_min} \\
        --codetection-exclude '${params.codetection_exclude}' \\
        --guardband-grid '${params.guardband_grid}' ${restrict} ${ipc} \\
        --ipc-target '${params.ipc_target}' \\
        --ipc-cq-max ${params.ipc_cq_max} \\
        --out-calls target_calls.csv \\
        --out-mosquito mosquito_level.csv \\
        --out-sensitivity crossamp_sensitivity.csv \\
        --out-ipc-summary ipc_summary.csv
    """
}

process REPLICATE_CONCORDANCE {
    input:
    path calls
    path assignment
    path round1_calls

    output:
    path 'replicate_pairs.csv',   emit: pairs
    path 'replicate_summary.csv', emit: summary

    script:
    def r1 = round1_calls.name != 'NO_ROUND1'
             ? "--round1-calls ${round1_calls}" : ''
    """
    replicate_concordance.py \\
        --calls ${calls} --assignment ${assignment} \\
        --targets '${params.targets}' ${r1} \\
        --out-pairs replicate_pairs.csv \\
        --out-summary replicate_summary.csv
    """
}

process RANDOMISATION_CHECK {
    input:
    path mosquito
    path wells
    path assignment
    path specimens

    output:
    path 'batch_diagnostic.csv',    emit: diagnostic
    path 'plate_positivity.csv',    emit: plate_positivity
    path 'allocation_balance.csv',  emit: balance

    script:
    """
    randomisation_check.py \\
        --mosquito ${mosquito} --wells ${wells} \\
        --assignment ${assignment} --specimens ${specimens} \\
        --strata '${params.balance_strata}' \\
        --icc-groups '${params.icc_groups}' \\
        --max-level-share ${params.max_level_share} \\
        --n-perm ${params.n_perm} --seed ${params.seed} \\
        --out-diagnostic batch_diagnostic.csv \\
        --out-plate plate_positivity.csv \\
        --out-balance allocation_balance.csv
    """
}

process REPORT {
    input:
    path mosquito
    path sensitivity
    path target_qc
    path control_qc
    path diagnostic
    path allocation_balance
    path ingest_audit
    path plate_map_report
    path replicate_summary
    path ipc_summary
    path cross_reactivity
    path reconciliation

    output:
    path 'analysis_report.md'
    path 'prevalence_by_zone.csv'

    script:
    """
    make_report.py --mosquito ${mosquito} --sensitivity ${sensitivity} \\
        --target-qc ${target_qc} --control-qc ${control_qc} \\
        --diagnostic ${diagnostic} \\
        --allocation-balance ${allocation_balance} \\
        --ingest-audit ${ingest_audit} \\
        --plate-map-report ${plate_map_report} \\
        --replicate-summary ${replicate_summary} \\
        --ipc-summary ${ipc_summary} \\
        --cross-reactivity ${cross_reactivity} \\
        --reconciliation ${reconciliation} \\
        --version 'v${workflow.manifest.version}' \\
        --out-report analysis_report.md \\
        --out-prevalence prevalence_by_zone.csv
    """
}
