"""
For a given set of WGAC coordinates, get the unique set of alignments for each
pair of sequences.
"""
import csv
import optparse
import sys


def get_unique_wgac(wgac_alignments_filename):
    unique_alignments = set()

    with open(wgac_alignments_filename, "r") as fh:
        wgac_alignments = csv.reader(fh, delimiter="\t")
        for row in wgac_alignments:
            alignment = tuple(row)

            # Reverse the current alignment by switching the coordinates of the
            # sequence pairs in the list while keeping the strand the same.
            reverse_alignment = tuple(alignment[4:-1] + alignment[3:4] + alignment[:3] + alignment[7:8])
            if reverse_alignment not in unique_alignments:
                unique_alignments.add(alignment)

    for alignment in unique_alignments:
        print "\t".join(alignment)


if __name__ == "__main__":
    parser = optparse.OptionParser(usage="%prog <wgac alignments>")
    options, args = parser.parse_args()

    if len(args) != 1:
        parser.print_usage()
        sys.exit(1)

    get_unique_wgac(args[0])
