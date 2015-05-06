#!/bin/env python
"""
Perform global alignment between all WGAC pairs and find single base mismatches
between each pair.

Approach:

1. Create a sequence file containing all sequences involved in WGAC pairwise
alignments. The first sequence in any alignment with "+" in the strand column
will be assigned the standard genomic orientation while the second sequence will
be reverse complemented.

2. Load sequence file with SeqIO index.

3. Load WGAC alignments.

4. For each WGAC alignment:
  a. Find sequences in index.
  b. Write sequences out to temporary files.
  c. Align sequences.
  d. Parse alignment for mismatches.
  e. Update mismatch coordinates for starting coordinate of each sequence accounting for reverse complement.
  f. Print genomic mismatch coordinates for both pairs.
"""
import argparse
import csv
from Bio import AlignIO, SeqIO
from Bio.Alphabet import generic_nucleotide
from Bio.Emboss.Applications import StretcherCommandline
from Bio.SeqRecord import SeqRecord
import logging
from mpi4py import MPI
import os
import sha
import shutil
import subprocess
import sys

# Setup logging.
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Log to stdout.
sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)
logger.addHandler(sh)

# Define MPI constants.
MASTER=0
TAG_WORK=11
TAG_DONE=12
TAG_DIE=13

# Default gap penalties from EMBOSS's stretcher alignment tool.
GAP_OPEN_COST=16
GAP_EXTEND_COST=4


def parse_wgac_alignment(wgac_alignment):
    """
    Parses a WGAC alignment from a list into a tuple of tuples for each sequence
    in the alignment.
    """
    wgac_alignment[1:3] = map(int, wgac_alignment[1:3])
    wgac_alignment[5:7] = map(int, wgac_alignment[5:7])
    return (
        tuple(wgac_alignment[:4]),
        tuple(wgac_alignment[4:])
    )


def get_alignment_name(alignment):
    alignment_name = "_".join(map(str, alignment))
    alignment_name = alignment_name.replace("+", "forward").replace("-", "reverse")
    return alignment_name


def align_sequences(asequence, bsequence, tmp_dir="/tmp"):
    """
    Aligns two sequences defined in the given filenames using EMBOSS's stretcher
    alignment tool.

    Returns a BioPython alignment object.
    """
    logger.debug("Align sequences %s and %s", asequence, bsequence)
    command = StretcherCommandline()
    command.asequence = asequence
    command.bsequence = bsequence
    command.outfile = ".".join(map(os.path.basename, (asequence, bsequence)) + ["stretcher"]).replace(".fasta", "")
    command.outfile = os.path.join(tmp_dir, command.outfile)
    command.aformat = "srspair"
    command.gapopen = GAP_OPEN_COST
    command.gapextend = GAP_EXTEND_COST
    stdout, stderr = command()

    logger.debug("Parse alignment for %s", command.outfile)
    alignment = AlignIO.read(command.outfile, "emboss")
    logger.debug("Loaded %i sequences from alignment", len(alignment))

    logger.debug("Removing alignment file %s", command.outfile)
    os.remove(command.outfile)

    return alignment


def convert_coordinates(coordinates, position):
    """
    Reverse complement a zero-based position coordinate if needed for the given
    range of coordinates including the strand of the coordinates.
    """
    if coordinates[3] == "-":
        length = coordinates[2] - coordinates[1]
        position = length - position - 1

    return position


# TODO: rename wgac_alignment to wgac_coordinates?
def parse_alignment_for_mismatches(wgac_alignment, alignment):
    """
    For a given WGAC alignment pair (coordinates) and a BioPython global
    alignment object of the two sequences in the WGAC pair, find all single base
    pair differences between the two sequences and return the genomic
    coordinates of the differences along with the bases that differ.
    """
    a_coordinates, b_coordinates = wgac_alignment
    a_sequence, b_sequence = alignment
    logger.debug("Sequence one: %s %s", a_sequence, a_coordinates)
    logger.debug("Sequence two: %s %s", b_sequence, b_coordinates)

    # Parse all differences between the two sequences (mismatches and indels).
    m = 0
    n = 0
    mismatches = []

    for i in xrange(alignment.get_alignment_length()):
        # If sequences aren't the same, there is a mismatch or gap.
        if a_sequence[i].upper() != b_sequence[i].upper():
            # Add the base for each sequence to the list of mismatches with
            # genomic coordinates. Take into account whether the sequence was
            # reverse complemented or not.
            a_start = a_coordinates[1] + convert_coordinates(a_coordinates, m)
            b_start = b_coordinates[1] + convert_coordinates(b_coordinates, n)

            # Track each mismatch pair in the score field of the BED output to
            # enable downstream processing. For instance, if someone wanted to
            # recalculate Ti/Tv ratio without double-counting they could do so
            # by mismatch number.
            mismatch_id = sha.sha("".join(map(str, (a_coordinates[0], a_start, b_coordinates[0], b_start)))).hexdigest()

            mismatches.append([a_coordinates[0], a_start, a_start + 1, ">".join((a_sequence[i], b_sequence[i])), mismatch_id, a_coordinates[3]])
            mismatches.append([b_coordinates[0], b_start, b_start + 1, ">".join((b_sequence[i], a_sequence[i])), mismatch_id, b_coordinates[3]])

        # Keep track of the current position in both sequences by excluding gaps
        # from the position count.
        if a_sequence[i] != "-":
            m += 1

        if b_sequence[i] != "-":
            n += 1

    logger.debug("Found %i mismatches", len(mismatches))
    return mismatches


def process_wgac_alignment(row, wgac_sequences, tmp_dir):
    """
    Process a given WGAC alignment row by extracting sequences in the WGAC pair,
    aligning sequences with stretcher, and parsing alignment for mismatches.
    """
    # Parse WGAC alignment.
    wgac_alignment = parse_wgac_alignment(row)
    logger.debug("Processing WGAC alignment: %s", wgac_alignment)

    # Create a temporary directory for this analysis.
    alignment_name = get_alignment_name(row)
    tmp_dir = os.path.join(tmp_dir, alignment_name)
    os.mkdir(tmp_dir)

    # Find sequences for alignment.
    sequences = [wgac_sequences.get("%s:%s-%s(%s)" % sequence)
                 for sequence in wgac_alignment]
    logger.debug("Found sequences: %s", sequences)

    # Write sequences to temporary files.
    sequence_filenames = []
    for sequence in sequences:
        base_sequence_name = sequence.id.replace(":", "_").replace("(-)", "_reverse").replace("(+)", "_forward").replace("-", "_")
        sequence_filenames.append(os.path.join(tmp_dir, "%s.fasta" % base_sequence_name))
        SeqIO.write([sequence], sequence_filenames[-1], "fasta")
        logger.debug("Wrote sequence to file %s", sequence_filenames[-1])

    # Align sequences.
    alignment = align_sequences(*sequence_filenames, tmp_dir=tmp_dir)

    # Parse alignment for mismatches.
    mismatches = parse_alignment_for_mismatches(wgac_alignment, alignment)

    # Delete temporary files.
    logger.debug("Delete temporary files in %s", tmp_dir)
    shutil.rmtree(tmp_dir)

    return mismatches


def run_locally(wgac_alignments_filename, wgac_sequences_filename, output_dir):
    """
    Locally runs global alignments of WGAC pairs and produces a BED file for
    mismatches in each alignment.
    """
    # Setup WGAC alignments index.
    tmp_dir = "/tmp/wgac"
    wgac_sequences = SeqIO.index(wgac_sequences_filename, "fasta")
    total_wgac_sequences = len(wgac_sequences)
    logger.info("Loaded %i WGAC sequences", total_wgac_sequences)

    # Load WGAC alignments into memory on the master node.
    with open(wgac_alignments_filename, "r") as fh:
        wgac_alignments = [row for row in csv.reader(fh, delimiter="\t")]
        logger.info("Loaded %i WGAC alignments", len(wgac_alignments))

    output_filename = "suns.bed"
    tmp_output = os.path.join(tmp_dir, output_filename)
    final_output = os.path.join(output_dir, output_filename)
    logger.debug("Write mismatches to tmp dir: %s", tmp_output)
    with open(tmp_output, "w") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator='\n')

        for alignment in wgac_alignments:
            mismatches = process_wgac_alignment(alignment, wgac_sequences, tmp_dir)
            writer.writerows(mismatches)

    # Sort BED files in place by chromosome and start position.
    return_value = subprocess.call("sort -k 1,1 -k 2,2n -u -o %s %s" % (tmp_output, tmp_output), shell=True)
    if return_value != 0:
        logger.error(
            "Sorting file %s on rank %s (node %s) returned non-zero exit code: %s",
            tmp_output,
            rank,
            node_name,
            return_value
        )

    logger.debug("Move mismatches to output dir: %s", final_output)
    shutil.move(tmp_output, final_output)
    logger.info("Done calculating SUNs")


def main(wgac_alignments_filename, wgac_sequences_filename, output_dir):
    """
    Distributes global alignments of WGAC pairs and produces a BED file for
    mismatches in each alignment.
    """
    # Setup MPI environment.
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    node_name = MPI.Get_processor_name()
    number_of_nodes = comm.Get_size()
    nodes = xrange(1, number_of_nodes)
    status = MPI.Status()

    # Setup WGAC alignments index.
    tmp_dir = "/tmp/wgac"
    wgac_sequences = SeqIO.index(wgac_sequences_filename, "fasta")
    total_wgac_sequences = len(wgac_sequences)
    logger.info("Loaded %i WGAC sequences", total_wgac_sequences)

    # Perform master/worker tasks.
    if rank == MASTER:
        # Send each alignment to a worker node.
        logger.info("Master node is %s", node_name)
        logger.info("Got %s nodes", number_of_nodes)

        # Load WGAC alignments into memory on the master node.
        with open(wgac_alignments_filename, "r") as fh:
            wgac_alignments = [row for row in csv.reader(fh, delimiter="\t")]
            logger.info("Loaded %i WGAC alignments", len(wgac_alignments))

        # Send work to nodes.
        working_nodes = []
        for node in nodes:
            if len(wgac_alignments) > 0:
                # Send work to each node as long as there is any left.
                work = wgac_alignments.pop(0)
                comm.send(work, dest=node, tag=TAG_WORK)
                working_nodes.append(node)
            else:
                # Kill any remaining nodes when there isn't any work left.
                comm.send("Die", dest=node, tag=TAG_DIE)

        # Send remaining work to nodes as they finish.
        while wgac_alignments:
            # Wait for any node to respond.
            result = comm.recv(
                source=MPI.ANY_SOURCE,
                tag=MPI.ANY_TAG,
                status=status
            )

            # Send next input to the node that finished.
            work = wgac_alignments.pop(0)
            comm.send(work, dest=status.source, tag=TAG_WORK)

        # Wait for working nodes to finish now that all work has been sent.
        for node in working_nodes:
            result = comm.recv(
                source=MPI.ANY_SOURCE,
                tag=MPI.ANY_TAG,
                status=status
            )
            logger.info("Node %s is done", status.source)

        # Shutdown all working nodes.
        for node in working_nodes:
            comm.send(None, dest=node, tag=TAG_DIE)

        logger.info("Master node is done.")
    else:
        # Non-master nodes perform the alignments until there is no more work.
        while True:
            alignment = comm.recv(
                source=MASTER,
                tag=MPI.ANY_TAG,
                status=status
            )

            if status.tag == TAG_DIE:
                logger.info("Rank %s (node %s) complete", rank, node_name)
                break
            else:
                logger.info("Rank %s (node %s) processing alignment: %s", rank, node_name, alignment)
                try:
                    mismatches = process_wgac_alignment(alignment, wgac_sequences, tmp_dir)

                    output_filename = "%s.bed" % get_alignment_name(alignment)
                    tmp_output = os.path.join(tmp_dir, output_filename)
                    final_output = os.path.join(output_dir, output_filename)
                    logger.debug("Write mismatches to tmp dir: %s", tmp_output)
                    with open(tmp_output, "w") as fh:
                        writer = csv.writer(fh, delimiter="\t", lineterminator='\n')
                        writer.writerows(mismatches)

                    # Sort BED files in place by chromosome and start position.
                    return_value = subprocess.call("sort -k 1,1 -k 2,2n -u -o %s %s" % (tmp_output, tmp_output), shell=True)
                    if return_value != 0:
                        logger.error(
                            "Sorting file %s on rank %s (node %s) returned non-zero exit code: %s",
                            tmp_output,
                            rank,
                            node_name,
                            return_value
                        )

                    logger.debug("Move mismatches to output dir: %s", final_output)
                    shutil.move(tmp_output, final_output)
                finally:
                    comm.send("Done", dest=MASTER, tag=TAG_DONE)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("wgac_alignments", help="tab-delimited file with pairwise alignments including relative orientation")
    parser.add_argument("wgac_sequences", help="FASTA file with one entry per distinct sequence in the pairwise alignments file"),
    parser.add_argument("output_dir", help="path where SUNs output files will be written")
    parser.add_argument("--local", action="store_true", help="calculate SUNs locally without distributing work through MPI")
    parser.add_argument("--verbose", action="store_true", help="display debugging log messages")
    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)
        sh.setLevel(logging.DEBUG)

    if not args.local:
        logger.info("Distributing SUNs calculation through MPI")
        main(args.wgac_alignments, args.wgac_sequences, args.output_dir)
    else:
        logger.info("Running SUNs calculation locally")
        run_locally(args.wgac_alignments, args.wgac_sequences, args.output_dir)
