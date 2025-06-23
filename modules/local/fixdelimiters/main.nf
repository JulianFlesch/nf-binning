process FIXDELIMITERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'biocontainers/ubuntu:22.04' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*cleaned.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outname = "${prefix}.cleaned.bed"
    """
    # Clean the input BED file to ensure consistent tab delimiters
    sed \\
        $args \\
        -E 's/\s+/\t/g' \\
        ${bed} \\
        > ${outname}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fixdelimiters: \$(sed --version | head -n 1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cleaned.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fixdelimiters: \$(sed --version | head -n 1)
    END_VERSIONS
    """
}
