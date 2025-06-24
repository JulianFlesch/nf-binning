#!/usr/bin/env python3

import os
import argparse
import csv


def parse_args():

    parser = argparse.ArgumentParser(description="Multiply the bedrgaph value and the overlap value from a bedfile")
    parser.add_argument("--version", action="version", version="%(prog)s 0.0.1")
    parser.add_argument("--input", help="Intersect file in bed format", type=str, required=True)
    parser.add_argument("--out", help="Ouput bed file location", default="normalized_intersect.bed")
    parser.add_argument("--overlap_col", type=int, default=-1, help="Index (0-based) describing in which column of the intersect bed file the overlap size is stored")
    parser.add_argument("--bedgraph_col", type=int, default=-2, help="Index (0-based) describing in which column (after running bedtools intersect) the bedgraph value is stored")

    args = parser.parse_args()

    return (args.input, args.out, args.overlap_col, args.bedgraph_col)


def mult_inter_bedgraph(bed_file: str, outfile: str, overlap_col: int, bedgraph_col: int):

    assert (os.path.exists(bed_file))
    assert (not os.path.extists(outfile))

    with open(bed_file, "r") as ib, open(outfile, "w") as ob:
        reader = csv.reaer(ib, delimiter="\t")
        writer = csv.writer(ob, delimiter="\t")

        for i, line in enumerate(reader):
            try:

                if i > 0 and len(line) == 0:
                    continue

                bedgraph_val = int(line[bedgraph_col])
                overlap_val = int(line[overlap_col])
                mult = bedgraph_val * overlap_val

                # write row with new normalized value
                line += [mult]
                writer.writerow(line)

            except IndexError:
                error_msg = "Invalid input file or overlap_column specified. Regions or bedgraph/overlap value could not be read! " + \
                            "(Input: %s, Linenum: %i, overlap_col: %i, bedgraph_col: %i)" % (bed_file, i, overlap_col, bedgraph_col)
                raise ValueError(error_msg)

            except ValueError:
                error_msg = "Possibly malformed input file: Unable to parse pos_start, pos_end, overlap or bedgraph value as integers. " + \
                            "(Input: %s, Linenum: %i)"
                raise ValueError(error_msg)


if __name__ == "__main__":
    input, output, overlap_col, bedgraph_col = parse_args()
    mult_inter_bedgraph(input, output, overlap_col, bedgraph_col)
