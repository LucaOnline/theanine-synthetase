"""Entrypoint analysis script."""

import json
from typing import List, Tuple

from dnds import dnds
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from alignment import AlignmentResult, align_sequences
from data_index import CSTSI, CSGSI, CSTSI_PROTEIN, CSGSI_PROTEIN
from dir_utils import get_data, make_output_dir, get_output
from dnds_prep import trim_for_dnds
from options import CLUSTER_COUNTS, DNDS_WINDOW_SIZES
from parse_fasta import parse_fasta
from sliding_window import sliding_window


def make_cluster_graphs(seq_filename: str, alignment_result: AlignmentResult):
    """
    Produces graphs with the mismatch cluster sizes and saves them to the
    output directory.
    """
    sns.set_theme()

    for clusters in CLUSTER_COUNTS:
        # Build DataFrame of mismatch windows
        df = pd.DataFrame(
            {"clusters": alignment_result.clustered_mismatches(cluster_count=clusters)}
        )

        sns.relplot(data=df["clusters"], kind="line")

        plt.title(f"mismatches per alignment cluster ({clusters} clusters)")
        plt.xlabel("cluster #")
        plt.ylabel("mismatches")

        fig = plt.gcf()
        fig.set_size_inches(7, 8)

        y_max = df["clusters"].max()
        plt.xticks(np.arange(clusters), np.arange(1, clusters + 1))
        plt.yticks(np.arange(y_max + 1), np.arange(y_max + 1))

        plt.savefig(get_output(f"{seq_filename}_clustered_mismatches_{clusters}.png"))


def make_dnds_graph(
    seq_filename: str,
    window_sizes: List[int],
    dnds_ratio_data: List[Tuple[List[int], List[float]]],
):
    """
    Produces several graphs showing dN/dS ratios across a whole sequence and saves
    them to the output directory. `window_sizes` and `dnds_ratio_data` are expected
    to be in the same order with respect to the analyses they represent.
    """
    sns.set_theme()

    for i, window_size in enumerate(window_sizes):
        df = pd.DataFrame({"bp": dnds_ratio_data[i][0], "ratio": dnds_ratio_data[i][1]})
        df = df.set_index("bp")

        sns.relplot(data=df, kind="line")

        plt.title(f"dN/dS ratios over windows of size {window_size}")
        plt.xlabel("window starting bp")
        plt.ylabel("dN/dS ratio")

        fig = plt.gcf()
        fig.set_size_inches(18, 6)

        plt.axhline(y=1.0, color="r", linestyle="--")

        plt.savefig(get_output(f"{seq_filename}_dnds_{window_size}.png"))


def sliding_window_dnds(
    sequence_1: str, sequence_2: str, window_size: int
) -> Tuple[List[int], List[float]]:
    """
    Performs a sliding-window dN/dS analysis over the provided sequences.
    Returns a list of tuples (start_base_pair, dnds_ratio).
    """
    windows = zip(
        sliding_window(sequence_1, n=window_size),
        sliding_window(sequence_2, n=window_size),
    )
    return (
        [i for i in range(len(sequence_1) - window_size + 1)],
        [dnds(*window_pair) for window_pair in windows],
    )


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

    if nucleotides:
        trimmed_alignment_1, trimmed_alignment_2 = trim_for_dnds(alignment_result)
        dnds_ratio_data = [
            sliding_window_dnds(trimmed_alignment_1, trimmed_alignment_2, window_size=i)
            for i in DNDS_WINDOW_SIZES
        ]
        make_dnds_graph(
            f"{seq2_filename}_dnds_ratios.png", DNDS_WINDOW_SIZES, dnds_ratio_data
        )

    for clusters in CLUSTER_COUNTS:
        clustered_mismatches = alignment_result.clustered_mismatches(
            cluster_count=clusters
        )
        clustered_mismatch_variance = alignment_result.clustered_mismatch_variance(
            cluster_count=clusters
        )

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
            )

        with open(get_output(f"{seq2_filename}_{clusters}.meta.json"), "w+") as f:
            json_output = {
                "percent_similarity": percent_similarity,
                "largest_mismatch_pos": largest_mismatch_pos,
                "largest_mismatch": largest_mismatch,
                "clustered_mismatch_variance": clustered_mismatch_variance,
                "clustered_mismatches": clustered_mismatches,
            }

            json.dump(
                json_output,
                f,
            )

    make_cluster_graphs(seq2_filename, alignment_result)


if __name__ == "__main__":
    analyze(CSTSI, CSGSI, nucleotides=True)
    analyze(CSTSI_PROTEIN, CSGSI_PROTEIN)
