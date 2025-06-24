process NORMALIZEOVERLAP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.norm_inter.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.norm_inter.bed"
    """
    normalize_overlap.py \\
        $args \\
        --input $bed \\
        --out $output \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        normalizeoverlap: \$(normalize_overlap.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.norm_inter.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        normalizeoverlap: \$(normalize_overlap.py --version)
    END_VERSIONS
    """
}
