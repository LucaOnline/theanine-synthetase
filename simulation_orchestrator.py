"""Entrypoint simulation orchestration script."""

import json
from subprocess import Popen, DEVNULL

from dir_utils import get_output
from monte_carlo import MonteCarloSimulationResult
from options import SIMULATION_COUNT


def orchestrate_simulations(parallelism: int = 1):
    """
    Divides the simulations to be done into `parallelism` chunks
    and runs them as separate Python processes. In general, the
    processes will be allocated uniformly to machine cores. This
    then aggregates the chunked simulation results into one final
    results file, `monte_carlo.agg.txt`. The `multiprocessing`
    `Pool` was not used here because the work to be allocated was
    too complicated to split into parts it can pickle.
    """

    if SIMULATION_COUNT % parallelism != 0:
        print("Simulations cannot be divided evenly into requested instance count!")
        exit(1)

    simulations = [
        Popen(
            [
                "python",
                "simulation.py",
                "--id",
                str(i),
                "--trials",
                str(SIMULATION_COUNT // parallelism),
            ],
            stdout=DEVNULL,
        )
        for i in range(parallelism)
    ]

    n_trials = 0
    n_successes = 0
    for i, simulation in enumerate(simulations):
        simulation.wait()
        with open(get_output(f"monte_carlo.{i}.json"), "r+") as f:
            data = json.load(f)
            n_trials += data["n_trials"]
            n_successes += data["n_successes"]

    p_value = n_successes / n_trials

    aggregated_result = MonteCarloSimulationResult(p_value, n_trials, n_successes)
    with open(get_output(f"monte_carlo.agg.txt"), "w+") as f:
        f.write(aggregated_result.format_result())


if __name__ == "__main__":
    # multiprocessing.cpu_count() can also be used, but there
    # are no guarantees that the CPU count will evenly divide
    # the number of trials.
    orchestrate_simulations(parallelism=8)
