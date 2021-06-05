import os
from random import shuffle
from typing import Callable

from alignment import align_sequences, AlignmentResult
from data_index import CSTSI, CSGSI
from monte_carlo import monte_carlo
from parse_fasta import parse_fasta


def get_data(filename: str) -> str:
    """Gets the relative file path for the provided file."""
    return "./data/" + filename


def make_output_dir():
    """Creates the output file directory, if it doesn't exist."""
    if not os.path.exists("./output"):
        os.mkdir("./output")


def get_output(filename: str) -> str:
    """Gets the relative file path for the provided output file."""
    return "./output/" + filename


def get_clustering_simulation_fn(
    cstsi_sequence: str, csgsi_sequence: str
) -> Callable[[], AlignmentResult]:
    """
    Creates a simulation function to be used in a Monte-Carlo simulation
    that takes a CsTSI and a CsGSI sequence and compares a shuffled
    CsGSI with CsTSI each time it is called.
    """

    def simulation_fn() -> AlignmentResult:
        rand_seq_list = list(csgsi_sequence)
        shuffle(rand_seq_list)
        rand_seq = "".join(rand_seq_list)

        return align_sequences(rand_seq, cstsi_sequence, nucleotides=True)

    return simulation_fn


def get_effect_size_fn(chunk_count: int) -> Callable[[AlignmentResult], float]:
    """
    Creates an effect size function for a Monte-Carlo simulation.
    The returned function returns the variance between the `chunk_count`
    chunks of the alignment result provided.
    """

    def effect_size_fn(alignment_result: AlignmentResult) -> float:
        return alignment_result.clustered_mismatch_variance(chunk_count)

    return effect_size_fn


if __name__ == "__main__":
    cstsi_mrna = list(parse_fasta(get_data(CSTSI)))[0][1]
    cluster_count = (
        15  # Number of clusters to use in mismatch clustering variance analysis
    )

    # Analyze CsGSI sequence
    print("Analyzing %s..." % CSGSI)

    csgsi_seq = list(parse_fasta(get_data(CSGSI)))[0][1]
    alignment_result = align_sequences(csgsi_seq, cstsi_mrna, nucleotides=True)

    print(
        "Variance between clusters: "
        + str(alignment_result.clustered_mismatch_variance(cluster_count=cluster_count))
    )

    largest_mismatch_pos, largest_mismatch = alignment_result.largest_mismatch()
    percent_similarity = 1 - (
        alignment_result.hamming_distance() / alignment_result.get_alignment_length()
    )

    # Output Supplementary Data 4
    make_output_dir()
    with open(get_output(CSGSI + ".aln.txt"), "w+") as f:
        f.write(alignment_result.format_result(line_length=100))
    with open(get_output(CSGSI + ".meta.txt"), "w+") as f:
        f.write(
            f"Percent similarity: {percent_similarity}\nLargest mismatch location: {largest_mismatch_pos}\nLargest mismatch size: {largest_mismatch}bp\n"
        )

    # Simulate random sequences
    simulation_result = monte_carlo(
        get_clustering_simulation_fn(cstsi_mrna, csgsi_seq),
        get_effect_size_fn(chunk_count=cluster_count),
        observed_effect_size=alignment_result.clustered_mismatch_variance(
            cluster_count=cluster_count
        ),
        n_trials=1000,
    )

    with open(get_output("monte_carlo.txt"), "w+") as f:
        f.write(simulation_result.format_result())

    simulation_result.examine()
