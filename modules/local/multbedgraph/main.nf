process MULTBEDGRAPH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when || meta.is_bedgraph

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.mult.bed"
    """
    mult_intersect_and_bedgraph.py \\
        $args \\
        --input $bed \\
        --out $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multbedgraph: \$(mult_intersect_and_bedgraph.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mult.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multbedgraph: \$(mult_intersect_and_bedgraph.py --version)
    END_VERSIONS
    """
}
