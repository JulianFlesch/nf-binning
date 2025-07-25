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


    // collect read_count as coverage statistic from all bed files
    ch_read_counts = Channel.empty()
    ch_sorted
        .map { meta, bed ->
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

            tuple(meta.id, tuple(num_reads, meta, bed))
        }
        .set { ch_read_counts }

    // collect total fragment length statistic from all bed files
    ch_fragment_len = Channel.empty()
    ch_sorted
        .map { meta, bed ->
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

                tuple(meta.id, total_fragment_len)
        }
        .set { ch_fragment_len }

    ch_stats = ch_fragment_len
        .join(ch_read_counts)
        .map{ id, fragment_len, read_counts_tuple -> 
            def (num_reads, meta, bed) = read_counts_tuple
            meta.total_fragment_len = fragment_len
            meta.fragment_count = num_reads
            tuple(meta, bed)
        }
    ch_stats.view()

    // Filter by read count 
    ch_filtered_1 = Channel.empty()
    if (params.fltr_min_num_reads) {
        ch_stats
            .filter { meta, bed ->
                // filter samples by read / fragment count
                meta.fragment_count >= params.fltr_min_num_reads
            }
            .set { ch_filtered_1 }
    } else {
        ch_filtered_1 = ch_stats
    }
    

    // Filter by total fragment length
    ch_filtered_2 = Channel.empty()
    if (params.fltr_min_tot_fragment_len) {
        ch_filtered_1
            .filter { meta, bed ->
                meta.total_fragment_len >= params.fltr_min_tot_fragment_len
            }
            .set { ch_filtered_2 }
    } else {
        ch_filtered_2 = ch_filtered_1
    }


    emit:
    sorted_bed      = ch_filtered_2            // channel: [ val(meta), [ bam ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
