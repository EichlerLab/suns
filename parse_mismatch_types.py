#!/bin/env python
"""
For a given BED file of SUNs with the base change in the name column (e.g.,
"A>G") and the unique id of each mismatch in the score column, calculate the
total number of each mismatch type including transitions, transversions,
insertions, and deletions.
"""
from argparse import ArgumentParser
from collections import defaultdict
import csv
import pprint


def get_mismatch_type(a_base, b_base):
    """
    For two given nucleotide values (A, T, C, G, and -) representing a mismatch
    or indel, return the type of mismatch. For substitutions, return transition
    or transversion. Otherwise, return insertion or deletion.
    """
    type_string = "".join(sorted((a_base, b_base)))

    if type_string in ("AG", "CT"):
        return "transition"
    elif type_string in ("AC", "AT", "CG", "GT"):
        return "transversion"
    elif a_base == "-":
        return "insertion"
    elif b_base == "-":
        return "deletion"
    else:
        raise Exception("Unrecognized mismatch type: %s" % type_string)


def parse_mismatch_types(suns_filename):
    mismatch_types = defaultdict(int)
    counted_mismatches = set()

    with open(suns_filename, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if row[4] not in counted_mismatches:
                mismatch_types[get_mismatch_type(*(row[3].split(">")))] += 1
                counted_mismatches.add(row[4])

    return mismatch_types


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("suns_filename", help="SUNs in BED file format with base change in name column and unique mismatch id in score column")
    args = parser.parse_args()

    mismatch_types = parse_mismatch_types(args.suns_filename)
    pprint.pprint(mismatch_types)
