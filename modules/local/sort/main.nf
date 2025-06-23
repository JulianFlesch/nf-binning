process SORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'biocontainers/ubuntu:22.04' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.sorted.bed"), emit: sorted
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outfile = "${prefix}.sorted.bed"
    """
    sort \\
        $args \\
        -k1,1 -k2,2n \\
        --field-separator=\\t \\
        $bed \\
        > $outfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | head -n 1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | head -n 1)
    END_VERSIONS
    """
}
