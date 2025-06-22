process DROPCOLUMNS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'biocontainers/ubuntu:22.04' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Simplify the bed file to only the region columns
    # and add a column to sum over containing "1" for each region
    cat $bed | cut -f 1,2,3 | sed s/$/\t1/g" > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dropcolumns: 0.0.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dropcolumns: 0.0.1
    END_VERSIONS
    """
}
