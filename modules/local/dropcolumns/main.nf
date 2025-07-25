process DROPCOLUMNS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'biocontainers/ubuntu:22.04' }"

    input:
    tuple val(meta), path(bed)
    val add_aggregation_column

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def last_column_val = (!meta.is_bedgraph && add_aggregation_column) ? "1" : "\$NF"
    """
    # Simplify the bed file to only the columns 1-3 (describing the region)
    # and the last one (describing overlap, bedgrapah, or similar)
    awk '{print \$1 "\t" \$2 "\t" \$3 "\t" $last_column_val}' $bed > ${prefix}.cut.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cut.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1)
    END_VERSIONS
    """
}
