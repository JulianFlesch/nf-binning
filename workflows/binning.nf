/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_binning_pipeline'
include { BEDTOOLS_INTERSECT     } from '../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_MAKEWINDOWS as BEDTOOLS_MAKEWINDOWS_500   } from '../modules/nf-core/bedtools/makewindows/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BINNING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_window_file
    bin_fixed_500

    main:

    ch_versions = Channel.empty()

    // Note: only runs, when a window file is provided
    ch_samplesheet
        .combine(ch_window_file)
        // Change the order of arguments to the bedtools process:
        // Intersection is calucated relative to the regions file
        .map { meta, bed, regions -> return [meta, regions, bed] }
        .set { ch_intersect }

    BEDTOOLS_INTERSECT(ch_intersect, tuple([], []))
    ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)


    if (bin_fixed_500)  {
        // Bin the bedfiles by regular regions, if window sizes are provided
        BEDTOOLS_MAKEWINDOWS_500(ch_samplesheet)
        ch_versions.mix(BEDTOOLS_MAKEWINDOWS_500.out.versions)
    }

    // * Preparation of Tumor-Normal files:
    // We computed Tumor-Normal pairs from the CCRE region files and sorted them by largest differences in absolute value.
    // >>Script scripts/TvsN.csh, absval.pl

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'binning_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
