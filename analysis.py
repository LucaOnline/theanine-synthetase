"""Entrypoint analysis script."""

import json

from dnds import dnds
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import seaborn as sns

from alignment import AlignmentResult, align_sequences
from data_index import CSTSI, CSGSI, CSTSI_PROTEIN, CSGSI_PROTEIN
from dir_utils import get_data, make_output_dir, get_output
from dnds_prep import trim_indels
from options import CLUSTER_COUNTS
from parse_fasta import parse_fasta


def make_cluster_graphs(seq_filename: str, alignment_result: AlignmentResult):
    """Produces graphs with the mismatch cluster sizes and saves them to the output directory."""
    sns.set_theme()

    for clusters in CLUSTER_COUNTS:
        # Build DataFrame of mismatch windows
        df = pd.DataFrame(
            {
                f"{clusters}_clusters": alignment_result.clustered_mismatches(
                    cluster_count=clusters
                )
            }
        )

        sns.relplot(data=df[f"{clusters}_clusters"], kind="line")

        plt.title(f"Mismatches per alignment cluster ({clusters} clusters)")
        plt.xlabel("Cluster #")
        plt.ylabel("Mismatches")

        fig = plt.gcf()
        fig.set_size_inches(7, 8)

        plt.xticks(np.arange(clusters), np.arange(1, clusters + 1))

        plt.savefig(get_output(f"{seq_filename}_clustered_mismatches_{clusters}.png"))


def analyze(seq1_filename: str, seq2_filename: str, nucleotides: bool = False):
    """
    Performs an alignment-based analysis on the sequences in the provided FASTA files.
    `nucleotides` should be `True` if the two filenames refer to nucleotide sequences.
    """
    # Read first sequence
    seq_name, seq1 = list(parse_fasta(get_data(seq1_filename)))[0]

    # Analyze second sequence
    print("Analyzing %s..." % seq2_filename)

    seq2 = list(parse_fasta(get_data(seq2_filename)))[0][1]
    alignment_result = align_sequences(seq2, seq1, nucleotides=nucleotides)

    largest_mismatch_pos, largest_mismatch = alignment_result.largest_mismatch()
    percent_similarity = 1 - (
        alignment_result.hamming_distance() / alignment_result.get_alignment_length()
    )

    for clusters in CLUSTER_COUNTS:
        clustered_mismatches = alignment_result.clustered_mismatches(
            cluster_count=clusters
        )
        clustered_mismatch_variance = alignment_result.clustered_mismatch_variance(
            cluster_count=clusters
        )

        if nucleotides:
            trimmed_alignment_1, trimmed_alignment_2 = trim_indels(alignment_result)
            dnds_ratio = dnds(trimmed_alignment_1, trimmed_alignment_2)

        # Output Supplementary Data 4
        make_output_dir()
        with open(get_output(f"{seq2_filename}_{clusters}.aln.txt"), "w+") as f:
            f.write(alignment_result.format_result(line_length=100))

        with open(get_output(f"{seq2_filename}_{clusters}.meta.txt"), "w+") as f:
            f.write(
                "Formatted metadata -- not for programmatic use.\n\n"
                + f"Information for alignment with {seq_name}:\n\n"
                + f"Percent similarity: {percent_similarity}\n"
                + f"Largest mismatch location: {largest_mismatch_pos}\n"
                + f"Largest mismatch size: {largest_mismatch}bp\n"
                + f"Variance between clusters ({clusters} clusters): {clustered_mismatch_variance}\n"
                + f"Clustered mismatches: {clustered_mismatches}\n"
                + (f"dN/dS ratio: {dnds_ratio}\n" if nucleotides else "")
            )

        with open(get_output(f"{seq2_filename}_{clusters}.meta.json"), "w+") as f:
            json_output = {
                "percent_similarity": percent_similarity,
                "largest_mismatch_pos": largest_mismatch_pos,
                "largest_mismatch": largest_mismatch,
                "clustered_mismatch_variance": clustered_mismatch_variance,
                "clustered_mismatches": clustered_mismatches,
            }

            if nucleotides:
                json_output["dnds_ratio"] = dnds_ratio

            json.dump(
                json_output,
                f,
            )

    make_cluster_graphs(seq2_filename, alignment_result)


if __name__ == "__main__":
    analyze(CSTSI, CSGSI, nucleotides=True)
    analyze(CSTSI_PROTEIN, CSGSI_PROTEIN)
