/*
 * Standard-curve and control QC, and the specificity evidence the calling
 * rule rests on.
 */

process PLATE_QC {
    input:
    path wells

    output:
    path 'target_qc.csv',  emit: target_qc
    path 'control_qc.csv', emit: control_qc

    script:
    """
    plate_qc.py --wells ${wells} \\
        --r2-min ${params.miqe_r2_min} \\
        --eff-min ${params.miqe_eff_min} --eff-max ${params.miqe_eff_max} \\
        --slope-min ${params.miqe_slope_min} --slope-max ${params.miqe_slope_max} \\
        --ctl-amplification-max ${params.ctl_amplification_max} \\
        --pos-ctl-cq-tolerance ${params.pos_ctl_cq_tolerance} \\
        --out-target target_qc.csv --out-control control_qc.csv
    """
}

process CROSS_REACTIVITY {
    input:
    path wells

    output:
    path 'cross_reactivity.json', emit: verdict
    path 'cross_reactivity.csv',  emit: detail

    script:
    """
    cross_reactivity.py --wells ${wells} \\
        --targets '${params.targets}' \\
        --cq-positive-max ${params.cq_positive_max} \\
        --out-verdict cross_reactivity.json \\
        --out-detail cross_reactivity.csv
    """
}
