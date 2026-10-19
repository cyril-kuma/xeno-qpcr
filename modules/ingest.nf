/*
 * Reading the plate exports, and confirming the plates that ran are the
 * plates the design specified.
 */

process INGEST_PLATE {
    tag "$plate_id"

    input:
    tuple val(plate_id), path(plate)

    output:
    path "wells_${plate_id}.csv",  emit: wells
    path "audit_${plate_id}.json", emit: audit

    script:
    def drift = params.allow_ingest_drift ? '--allow-ingest-drift' : ''
    """
    ingest_plate.py \\
        --plate ${plate} \\
        --sheet '${params.sheet}' \\
        --header-row ${params.header_row} \\
        --validation-plate-pattern '${params.validation_plate_pattern}' \\
        --out wells_${plate_id}.csv \\
        --audit audit_${plate_id}.json ${drift}
    """
}

process MERGE_WELLS {
    input:
    path wells
    path audits

    output:
    path 'wells_merged.csv', emit: wells
    path 'ingest_audit.csv', emit: audit

    script:
    """
    merge_ingest.py --wells ${wells} --audits ${audits} \\
        --out-wells wells_merged.csv --out-audit ingest_audit.csv
    """
}

/*
 * The gate. Every stage after this one consumes `wells`, which this process
 * emits only once the specimen wells match the assignment, so the dataflow
 * itself enforces the order rather than a convention.
 */
process VALIDATE_PLATE_MAPS {
    input:
    path wells
    path plate_map

    output:
    path 'wells.csv',                 emit: wells
    path 'plate_map_validation.csv',  emit: report

    script:
    def allow = params.allow_plate_map_mismatch ? '--allow-plate-map-mismatch' : ''
    """
    validate_plate_maps.py \\
        --wells ${wells} --plate-map ${plate_map} \\
        --out-wells wells.csv \\
        --out-report plate_map_validation.csv ${allow}
    """
}
