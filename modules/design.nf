/*
 * Design-time stages: the allocation, and the sample sheet the analysis reads.
 */

process RANDOMISE {
    tag "seed ${params.randomisation_seed}"

    input:
    path specimens
    path first_round_slots

    output:
    path 'sample_assignment_new.csv', emit: assignment
    path 'allocation_balance.csv',    emit: balance

    script:
    def slots = first_round_slots.name != 'NO_FIRST_ROUND'
                ? "--first-round-slots ${first_round_slots}" : ''
    """
    randomise.py \\
        --specimens ${specimens} ${slots} \\
        --n-plates ${params.n_plates} \\
        --slots-per-plate ${params.slots_per_plate} \\
        --seed ${params.randomisation_seed} \\
        --out-assignment sample_assignment_new.csv \\
        --out-balance allocation_balance.csv
    """
}

process IMPORT_SAMPLE_SHEET {
    tag "${sample_sheet.name}"

    input:
    path sample_sheet

    output:
    path 'specimens.csv',                emit: specimens
    path 'sample_assignment.csv',        emit: assignment
    path 'plate_map_expected.csv',       emit: plate_map
    path 'specimen_reconciliation.csv',  emit: reconciliation

    script:
    """
    import_sample_sheet.py \\
        --sample-sheet '${sample_sheet}' \\
        --max-unmatched ${params.max_unmatched_specimens} \\
        --out-specimens specimens.csv \\
        --out-assignment sample_assignment.csv \\
        --out-plate-map plate_map_expected.csv \\
        --out-reconciliation specimen_reconciliation.csv
    """
}
