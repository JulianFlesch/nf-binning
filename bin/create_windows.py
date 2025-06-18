 #!/bin/python

import argparse
import csv
import copy
from dataclasses import dataclass
from typing import Self


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
        print("updated start:", self.start)

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


def merge_into_windows(input_bed: str, window_size: int) -> list[tuple[str, int, int]]:
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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Merge regions in a BED file to regions that are multiples of a specified length.")

    input_bed = "assets/testdata/foo_0.bed"
    window_size = 500
    output_bed = "windows.bed"

    windows = merge_into_windows(input_bed, window_size)

    for i, win in enumerate(windows):
        print(f"{i}\t{win[0]}\t{win[1]}\t{win[2]}")

        if i == 10:
            break

