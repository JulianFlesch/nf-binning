#!/usr/bin/env python3

import os
import argparse
import csv


def parse_args():

    parser = argparse.ArgumentParser(description="Normalize the overlap column in intersected bed files by the region length")
    parser.add_argument("--version", action="version", version="%(prog)s 0.0.1")
    parser.add_argument("--input", help="Intersect file in bed format", type=str, required=True)
    parser.add_argument("--out", help="Ouput bed file location", default="normalized_intersect.bed")
    parser.add_argument("--overlap_col", type=int, default=-1, help="Index (0-based) describing in which column of the intersect bed file the overlap size is stored")
    parser.add_argument("--round", type=int, default=3, help="Number of digits to round to.")

    args = parser.parse_args()

    return (args.input, args.out, args.overlap_col, args.round)


def normalize_overlap(bed_file: str, outfile: str, overlap_col: int, round_to: int=3):

    assert (os.path.exists(bed_file))
    assert (not os.path.exists(outfile))

    with open(bed_file, "r") as ib, open(outfile, "w") as ob:
        reader = csv.reader(ib, delimiter="\t")
        writer = csv.writer(ob, delimiter="\t")

        for i, line in enumerate(reader):
            try:

                if i > 0 and len(line) == 0:
                    continue

                pos_start = int(line[1])
                pos_end = int(line[2])
                overlap = int(line[overlap_col])
                normalized = round(
                                min(1, overlap / (pos_end - pos_start)),  # values should never be bigger than 1
                                round_to
                            )

                # write row with new normalized value
                line[overlap_col] = normalized
                writer.writerow(line)

            except IndexError:
                error_msg = "Invalid input file or overlap_column specified. Regions or overlap value could not be read! " + \
                            "(Input: %s, Linenum: %i, overlap_col: %i)" % (bed_file, i, overlap_col)
                raise ValueError(error_msg)

            except ValueError:
                error_msg = "Possibly malformed input file: Unable to parse pos_start, pos_end, or overlap as integers. " + \
                            "(Input: %s, Linenum: %i)"
                raise ValueError(error_msg)


if __name__ == "__main__":
    input, output, overlap_col, round_to = parse_args()
    normalize_overlap(input, output, overlap_col, round_to)
