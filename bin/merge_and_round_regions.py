#!/usr/bin/env python3

import argparse
import csv
import copy
from dataclasses import dataclass
from typing import Self
import os

@dataclass
class Window:
    chrom: str|None = None
    start: int|None = None
    end: int|None = None

    @property
    def is_empty(self) -> bool:
        return self.chrom is None and self.start is None and self.end is None

    @property
    def has_chrom(self) -> bool:
        return self.chrom is not None

    @property
    def has_start(self) -> bool:
        return self.start is not None

    @property
    def has_end(self) -> bool:
        return self.end is not None

    def set_start(self, start: int, round_to=1) -> None:
        # Update the value of start, rounding down to the nearest multiple of round_to
        self.start = start // round_to * round_to

    def set_end(self, end:int, round_to=1) -> None:
        # Update the value of end, rounding up to the nearest multiple of round_to
        self.end = (end // round_to * round_to) + round_to

    def update(self, window: Self, round_to=1) -> None:
        self.chrom = window.chrom
        self.set_start(window.start, round_to=round_to)
        self.set_end(window.end, round_to=round_to)

    def reset(self):
        self.chrom = None
        self.start = None
        self.end = None

    def to_tuple(self) -> tuple[str, int, int]:
        return (self.chrom, self.start, self.end)

    def __str__(self) -> str:
        return ", ".join([str(e) for e in self.to_tuple()])


def merge_and_round_regions(input_bed: str, window_size: int) -> list[tuple[str, int, int]]:
    """
    Rounds region dimensions by window_size and merges overlapping regions after rounding.
    """
    out_regions = []
    with open(input_bed, "r") as fb:
        reader = csv.reader(fb, delimiter='\t')
        current_window = Window()
        for line in reader:
            new_window = Window(line[0], int(line[1]), int(line[2]))

            if current_window.is_empty:
                current_window.update(new_window, round_to=window_size)

            elif current_window.chrom != new_window.chrom:
                out_regions.append(current_window.to_tuple())
                current_window.update(new_window, round_to=window_size)

            else:

                if not current_window.has_start:
                    current_window.set_start(new_window.start, round_to=window_size)

                # Save and start new region if the start is beyond the current end
                elif current_window.has_end and new_window.start > current_window.end + window_size:
                    out_regions.append(current_window.to_tuple())
                    current_window.update(new_window, round_to=window_size)

                # Extend the current region if the end is inside the current window
                else:
                    current_window.set_end(new_window.end, round_to=window_size)

    return out_regions


def parse_args():

    parser = argparse.ArgumentParser(description="Merge bed regions into multiples of fixed-sized windows sizes.")

    parser.add_argument("--input_bed", "-i", type=str, required=True, help="Input BED file with regions to create fixed-sized windows over.")
    parser.add_argument("--window_size", "-w", type=int, required=True, help="size of the windows to create, in base pairs.")
    parser.add_argument("--output_bed", "-o", type=str, default="windows.bed", help="Output BED file to write the windows to.")
    args = parser.parse_args()

    return args.input_bed, args.window_size, args.output_bed


if __name__ == "__main__":

    input_bed, window_size, output_bed = parse_args()

    # Check if the output file exists
    if os.path.exists(output_bed):
        raise ValueError("Output file already exists")

    # ensure the input BED file exists
    if not os.path.exists(input_bed):
        raise FileNotFoundError(f"Input BED file {input_bed} does not exist.")

    windows = merge_and_round_regions(input_bed, window_size)

    # Write the windows to the output BED file
    with open(output_bed, "w") as ob:
        writer = csv.writer(ob, delimiter='\t')
        for window in windows:
            writer.writerow(window)

    print(f"Created {len(windows)} regions from {input_bed} and saved to {output_bed}.")
