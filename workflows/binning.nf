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
include { BEDTOOLS_MERGE } from '../modules/nf-core/bedtools/merge/main'
include { SIMPLIFY_REGIONS } from '../modules/local/simplify_regions/main'
include { DROPCOLUMNS as DROPCOLUMNS_REGIONS } from '../modules/local/dropcolumns/main'
include { DROPCOLUMNS as DROPCOLUMNS_WINDOWS } from '../modules/local/dropcolumns/main'
include { BEDTOOLS_GROUPBY as BEDTOOLS_GROUPBY_REGIONS } from '../modules/nf-core/bedtools/groupby/main'
include { BEDTOOLS_GROUPBY as BEDTOOLS_GROUPBY_WINDOWS } from '../modules/nf-core/bedtools/groupby/main'
include { CAT_SORT } from '../modules/local/cat_sort/main'
include { NORMALIZEOVERLAP as NORMALIZEOVERLAP_REGIONS } from '../modules/local/normalizeoverlap/main'
include { NORMALIZEOVERLAP as NORMALIZEOVERLAP_WINDOWS } from '../modules/local/normalizeoverlap/main'
include { MULTBEDGRAPH as MULTBEDGRAPH_REGIONS } from '../modules/local/multbedgraph/main'
include { MULTBEDGRAPH as MULTBEDGRAPH_WINDOWS } from '../modules/local/multbedgraph/main'
include { SORTANDFILTER } from '../subworkflows/local/sortandfilter'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BINNING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_regions_file
    bin_fixed_500

    main:

    ch_versions = Channel.empty()

    // Specifies if in DROPCOLUMNS a column containing only "1" is added
    add_aggregation_col = !(params.use_bedgraph_value || params.normalize_overlap)

    SORTANDFILTER(ch_samplesheet)
    ch_versions.mix(SORTANDFILTER.out.versions)

    // BIN BY PREDEFINED REGIONS
    // -------------------------
    if (params.regions_file) {
        // (Note: only executed, when a window file is provided in ch_regions_file)
        SORTANDFILTER.out.sorted_bed
            .combine(ch_regions_file)
            // Change the order of arguments to the bedtools process:
            // Intersection is calucated relative to the regions file
            .map { meta, bed, regions -> return [meta, regions, bed] }
            .set { ch_intersect }

        BEDTOOLS_INTERSECT_REGIONS(ch_intersect, tuple([], []))
        ch_versions.mix(BEDTOOLS_INTERSECT_REGIONS.out.versions)
        ch_bedtools_intersect_regions = BEDTOOLS_INTERSECT_REGIONS.out.intersect

        // Multiply Bedgraph Value by Overlap
        // only runs for bedgraph files or if use_bedgraph_value is set
        MULTBEDGRAPH_REGIONS(ch_bedtools_intersect_regions)
        ch_versions.mix(MULTBEDGRAPH_REGIONS.out.versions)
        ch_mult_regions = MULTBEDGRAPH_REGIONS.out.bed

        // Drop all columns except 1,2,3 and the last (if either bedgraph value or noramlized overlap is given)
        // Add a column of for summing with groupby
        DROPCOLUMNS_REGIONS(ch_mult_regions, add_aggregation_col)
        ch_versions.mix(DROPCOLUMNS_REGIONS.out.versions)
        ch_dropcol_regions = DROPCOLUMNS_REGIONS.out.bed

        // Bedtools groupby and sum!
        BEDTOOLS_GROUPBY_REGIONS(ch_dropcol_regions, 4)
        ch_versions.mix(BEDTOOLS_GROUPBY_REGIONS.out.versions)
        ch_bedtools_groupby_regions = BEDTOOLS_GROUPBY_REGIONS.out.bed

        // Normalize and replace overlap value from bedtools intersect
        ch_normalize_inter_regions = Channel.empty()
        if (params.normalize_overlap) {
            NORMALIZEOVERLAP_REGIONS(ch_bedtools_groupby_regions)
            ch_versions.mix(NORMALIZEOVERLAP_REGIONS.out.versions)
            ch_normalize_inter_regions = NORMALIZEOVERLAP_REGIONS.out.bed
        } else {
            ch_normalize_inter_regions = ch_bedtools_groupby_regions
        }

    }

    // BIN BY FIXED 500bp REGIONS
    // --------------------------
    window_size = 500
    if (bin_fixed_500) {

        SORTANDFILTER.out.sorted_bed
            .map { _meta, bed -> return bed }
            .collect()
            .set { ch_beds }

        ch_dummy_meta = Channel.value([id: "all_samples", condition: 0])

        ch_dummy_meta
            .combine(ch_beds.toList())
            .set { ch_concat }

        CAT_SORT(ch_concat)
        ch_versions.mix(CAT_SORT.out.versions)

        // Merge overlapping regions
        BEDTOOLS_MERGE(CAT_SORT.out.sorted)
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

        // only for bedgraph files or if use_bedgraph_value is set
        MULTBEDGRAPH_WINDOWS(BEDTOOLS_INTERSECT_WINDOWS.out.intersect)
        ch_versions.mix(MULTBEDGRAPH_WINDOWS.out.versions)
        ch_mult_windows = MULTBEDGRAPH_WINDOWS.out.bed

        // Drop all columns except 1,2,3. Add a column containing "1" for each region
        DROPCOLUMNS_WINDOWS(ch_mult_windows, add_aggregation_col)
        ch_versions.mix(DROPCOLUMNS_WINDOWS.out.versions)
        ch_dropcol_windows = DROPCOLUMNS_WINDOWS.out.bed

        // bedtools groupby and sum
        BEDTOOLS_GROUPBY_WINDOWS(ch_dropcol_windows, 4)
        ch_versions.mix(BEDTOOLS_GROUPBY_WINDOWS.out.versions)
        ch_bedtools_groupby_windows = BEDTOOLS_GROUPBY_WINDOWS.out.bed

        // Calculated Overlap normalized by window size
        ch_normalize_windows = Channel.empty()
        if (params.normalize_overlap) {
            NORMALIZEOVERLAP_WINDOWS(ch_bedtools_groupby_windows)
            ch_versions.mix(NORMALIZEOVERLAP_WINDOWS.out.versions)
            ch_normalize_windows = NORMALIZEOVERLAP_WINDOWS.out.bed
        } else {
            ch_normalize_windows = ch_bedtools_groupby_windows
        }
    }

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
