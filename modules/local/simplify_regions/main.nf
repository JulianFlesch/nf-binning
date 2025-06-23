
process SIMPLIFY_REGIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python=3.11':
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(bed)
    val windows_size

    output:
    tuple val(meta), path("merged_and_rounded.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    merge_and_round_regions.py \\
        --input_bed ${bed} \\
        --window_size ${windows_size} \\
        --output_bed merged_and_rounded.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_and_round_regions: \$(merge_and_roud-regions.py --version )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch windows.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_and_round_regions: \$(merge_and_roud-regions.py --version )
    END_VERSIONS
    """
}
