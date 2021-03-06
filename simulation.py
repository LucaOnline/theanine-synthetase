"""Entrypoint simulation script."""

import argparse
from random import shuffle
from typing import Callable

from alignment import align_sequences, AlignmentResult
from data_index import CSTSI_PROTEIN, CSGSI_PROTEIN
from dir_utils import get_data, make_output_dir, get_output
from monte_carlo import monte_carlo
from options import CLUSTER_COUNTS
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

        return align_sequences(rand_seq, cstsi_sequence, nucleotides=False)

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
    cstsi_seq = list(parse_fasta(get_data(CSTSI_PROTEIN)))[0][1]

    # Analyze CsGSI sequence
    print("Analyzing %s..." % CSGSI_PROTEIN)

    csgsi_seq = list(parse_fasta(get_data(CSGSI_PROTEIN)))[0][1]
    alignment_result = align_sequences(csgsi_seq, cstsi_seq, nucleotides=False)

    for clusters in CLUSTER_COUNTS:
        print(
            f"Variance between clusters ({clusters} clusters): {str(alignment_result.clustered_mismatch_variance(cluster_count=clusters))}"
        )

        # Simulate random sequences
        simulation_result = monte_carlo(
            get_clustering_simulation_fn(cstsi_seq, csgsi_seq),
            get_effect_size_fn(cluster_count=clusters),
            observed_effect_size=alignment_result.clustered_mismatch_variance(
                cluster_count=clusters
            ),
            n_trials=n_trials,
        )

        make_output_dir()
        with open(
            get_output(f"monte_carlo_{clusters}.{simulation_id}.json"), "w+"
        ) as f:
            f.write(simulation_result.to_json())

        simulation_result.examine()
