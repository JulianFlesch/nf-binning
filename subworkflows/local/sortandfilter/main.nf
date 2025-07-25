include { SORT } from '../../../modules/local/sort/main'
include { FIXDELIMITERS } from '../../../modules/local/fixdelimiters/main'


workflow SORTANDFILTER {

    take:
    ch_bed // channel: [ val(meta), [ bam ] ]

    main:
    ch_versions = Channel.empty()

    ch_fixdelimiters = Channel.empty()
    if (params.fix_bedfile_delimiters) {
        FIXDELIMITERS(ch_bed)
        ch_versions.mix(FIXDELIMITERS.out.versions)
        ch_fixdelimiters = FIXDELIMITERS.out.bed
    } else {
        ch_fixdelimiters = ch_bed
    }

    ch_sorted = Channel.empty()
    if (params.sort_inputs) {
        SORT(ch_fixdelimiters)
        ch_versions.mix(SORT.out.versions)
        ch_sorted = SORT.out.sorted
    } else {
        ch_sorted = ch_fixdelimiters
    }

    // Filter by read count 
    ch_filtered_1 = Channel.empty()
    if (params.fltr_min_num_reads) {
        ch_sorted
            .filter { meta, bed ->
                def num_reads

                // Note: for bedgraphs with coverage values, this might overestimate the read count,
                // but it is the only actionable information.

                // TODO: should we always just count the lines?
                if (meta.is_bedgraph) {
                    // For bedGraph files, sum the coverage value expexted in the last column
                    // bed.readLines() reads the file into a list of strings (one per line)
                    num_reads = bed.readLines().sum { line ->
                        // Trim whitespace, split by any whitespace, and take the last element
                        line.trim().split(/\s+/)[-1] as long
                    }
                } else {
                    // For regular BED files, simply count the number of lines
                    num_reads = bed.readLines().size()
                }

                // filter samples by read / fragment count
                num_reads >= params.fltr_min_num_reads
            }
            .set { ch_filtered_1 }
    } else {
        ch_filtered_1 = ch_sorted
    }

    // Filter by total fragment length
    ch_filtered_2 = Channel.empty()
    if (params.fltr_min_tot_fragment_len) {
        ch_filtered_1
            .filter { meta, bed ->
                def total_fragment_len
                def splitted
                def start
                def end

                // iterate over the bedfile, summing the regions sizes
                // to obtain the total fragment length
                total_fragment_len = bed.readLines().sum { line ->
                    // Trim whitespace, split by any whitespace
                    splitted = line.trim().split(/\s+/)
                    end = splitted[2] as long
                    start = splitted[1] as long

                    // For bedgraph files that have a coverage column,
                    // multiply the region size
                    if (meta.is_bedgraph && splitted.size() > 3) {
                        (end - start) * splitted[-1]

                    // otherwise just calculate region size
                    } else {
                        end - start
                    }
                }

                // filter samples by total fragment length
                total_fragment_len >= params.fltr_min_tot_fragment_len
            }
            .set { ch_filtered_2 }
        ch_filtered_2.view()
    } else {
        ch_filtered_2 = ch_filtered_1
    }


    emit:
    sorted_bed      = ch_filtered_2            // channel: [ val(meta), [ bam ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
