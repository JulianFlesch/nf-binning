/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_binning_pipeline'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_REGIONS } from '../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_WINDOWS } from '../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_MAKEWINDOWS as BEDTOOLS_MAKEWINDOWS_500   } from '../modules/nf-core/bedtools/makewindows/main'
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
include { BEDTOOLS_SORT } from '../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_MERGE } from '../modules/nf-core/bedtools/merge/main'
include { SIMPLIFY_REGIONS } from '../modules/local/simplify_regions/main'
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


    // BIN BY PREDEFINED REGIONS
    // -------------------------
    // (Note: only executed, when a window file is provided in ch_window_file)
    ch_samplesheet
        .combine(ch_window_file)
        // Change the order of arguments to the bedtools process:
        // Intersection is calucated relative to the regions file
        .map { meta, bed, regions -> return [meta, regions, bed] }
        .set { ch_intersect }

    BEDTOOLS_INTERSECT_REGIONS(ch_intersect, tuple([], []))
    ch_versions.mix(BEDTOOLS_INTERSECT_REGIONS.out.versions)


    // BIN BY FIXED 500bp REGIONS
    // --------------------------
    window_size = 500
    if (bin_fixed_500) {

        ch_samplesheet
            .map { _meta, bed -> return bed }
            .collect()
            .set { ch_beds }

        ch_dummy_meta = Channel.value([id: "all_beds", condition: 0])

        ch_dummy_meta
            .combine(ch_beds.toList())
            .set { ch_concat }

        // Concat all bed files from the samplesheet
        CAT_CAT(ch_concat)
        ch_versions.mix(CAT_CAT.out.versions)

        // Sort the concatenated bed file
        BEDTOOLS_SORT(CAT_CAT.out.file_out, [])
        ch_versions.mix(BEDTOOLS_SORT.out.versions)

        // Merge overlapping regions
        BEDTOOLS_MERGE(BEDTOOLS_SORT.out.sorted)
        ch_versions.mix(BEDTOOLS_MERGE.out.versions)

        // Round and simplify the merged regions to the window size
        SIMPLIFY_REGIONS(BEDTOOLS_MERGE.out.bed, window_size)
        ch_versions.mix(SIMPLIFY_REGIONS.out.versions)

        // Create a bedfile with regular regions, if window sizes are provided
        BEDTOOLS_MAKEWINDOWS_500(SIMPLIFY_REGIONS.out.bed)
        ch_versions.mix(BEDTOOLS_MAKEWINDOWS_500.out.versions)

        // Intersect the window created window file with the bed files from our samplesheet
        BEDTOOLS_MAKEWINDOWS_500.out.bed
            .map { _meta, bed -> return bed }
            .combine(ch_samplesheet)
            // Change the order of arguments to the bedtools process:
            // Intersection is calucated relative to the regions file
            .map { regions, meta, bed -> return [meta, regions, bed] }
            .set { ch_intersect_windows }

        BEDTOOLS_INTERSECT_WINDOWS(ch_intersect_windows, tuple([], []))
        ch_versions.mix(BEDTOOLS_INTERSECT_WINDOWS.out.versions)

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
