import argparse
from random import shuffle
import sys
from typing import Callable

from alignment import align_sequences, AlignmentResult
from data_index import CSTSI, CSGSI
from dir_utils import get_data, make_output_dir, get_output
from monte_carlo import monte_carlo
from options import CLUSTER_COUNT
from parse_fasta import parse_fasta


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


def get_effect_size_fn(cluster_count: int) -> Callable[[AlignmentResult], float]:
    """
    Creates an effect size function for a Monte-Carlo simulation.
    The returned function returns the variance between the `cluster_count`
    chunks of the alignment result provided.
    """

    def effect_size_fn(alignment_result: AlignmentResult) -> float:
        return alignment_result.clustered_mismatch_variance(cluster_count)

    return effect_size_fn


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Run a Monte-Carlo simulation for CsTSI and CsGSI alignment."
    )
    parser.add_argument("--id", dest="simulation_id", type=int)
    parser.add_argument("--trials", dest="n_trials", type=int)
    args = parser.parse_args()

    simulation_id = args.simulation_id
    n_trials = args.n_trials

    # Read CsTSI sequence
    cstsi_mrna = list(parse_fasta(get_data(CSTSI)))[0][1]

    # Analyze CsGSI sequence
    print("Analyzing %s..." % CSGSI)

    csgsi_seq = list(parse_fasta(get_data(CSGSI)))[0][1]
    alignment_result = align_sequences(csgsi_seq, cstsi_mrna, nucleotides=True)

    print(
        "Variance between clusters: "
        + str(alignment_result.clustered_mismatch_variance(cluster_count=CLUSTER_COUNT))
    )

    # Simulate random sequences
    simulation_result = monte_carlo(
        get_clustering_simulation_fn(cstsi_mrna, csgsi_seq),
        get_effect_size_fn(cluster_count=CLUSTER_COUNT),
        observed_effect_size=alignment_result.clustered_mismatch_variance(
            cluster_count=CLUSTER_COUNT
        ),
        n_trials=n_trials,
        verbose=True,
    )

    make_output_dir()
    with open(get_output(f"monte_carlo.{simulation_id}.txt"), "w+") as f:
        f.write(simulation_result.format_result())

    simulation_result.examine()
