import json
from math import nan
from time import perf_counter
from typing import Callable, TypeVar

T = TypeVar("T")


class MonteCarloSimulationResult:
    """
    MonteCarloSimulationResult represents the result of a Monte-Carlo simulation.
    """

    def __init__(self, p_value: float, n_trials: int, n_successes: int):
        """
        Produces a new MonteCarloSimulationResult representing the result of a Monte-Carlo simulation.
        """
        self.p_value = p_value
        self.n_trials = n_trials
        self.n_successes = n_successes

    def examine(self):
        """Prints information about the simulation to the console."""
        print(self.format_result())

    def format_result(self):
        """Formats the results of the simulation for external use."""
        return f"Monte-Carlo simulation summary:\nTrials: {self.n_trials}\nSuccesses: {self.n_successes}\np-value: {self.p_value}"

    def to_json(self):
        """Formats the results of the simulation as a JSON string and returns it."""
        return json.dumps(
            {
                "p_value": self.p_value,
                "n_trials": self.n_trials,
                "n_successes": self.n_successes,
            }
        )

    def get_p_value(self):
        """Returns the estimated p-value of the simmulation."""
        return self.p_value

    def get_trial_count(self):
        """Returns the number of trials in the simmulation."""
        return self.n_trials

    def get_success_count(self):
        """Returns the number of successes in the simmulation."""
        return self.n_successes


def monte_carlo(
    simulation_fn: Callable[[], T],
    effect_size_fn: Callable[[T], float],
    observed_effect_size: float,
    n_trials: int,
    verbose: bool = False,
) -> MonteCarloSimulationResult:
    """
    Runs a Monte-Carlo simulation. Returns an object representing the result of the simulation,
    which will contain the final p-value.

    `simulation_fn`: A function that runs a single trial of the simulation.
    `effect_size_fn`: A function that calculates an effect size from the return value of the simulation function.
    `observed_effect_size`: The observed effect size that is being tested.
    `n_trials`: The number of trials to run. More trials will more precisely estimate the p-value.
    `verbose`: Whether or not to perform logging to the console while the simulation is running.
    """

    if verbose:
        print("Beginning Monte-Carlo simulation with %d trials." % n_trials)
        start_time = perf_counter()

    def get_next_result(trial: int) -> bool:
        if verbose:
            elapsed_time = perf_counter() - start_time
            eta = (elapsed_time / trial) * (n_trials + 1 - trial) if trial > 1 else nan
            eta_units = "seconds"

            if eta > 59:
                eta /= 60
                eta_units = "minutes"

            if eta > 59:
                eta /= 60
                eta_units = "hours"

            if eta > 23:
                eta /= 24
                eta_units = "days"

            print(
                "Running simulation %d/%d, estimated %.3f %s remaining..."
                % (trial, n_trials, eta, eta_units)
            )

        next_effect_size = effect_size_fn(simulation_fn())
        next_result = next_effect_size >= observed_effect_size

        if verbose:
            if next_result:
                print("Simulation %d suceeded." % trial)
            else:
                print("Simulation %d failed." % trial)

        return next_result

    trials = [get_next_result(trial + 1) for trial in range(n_trials)]
    n_successes = sum(trials)  # Yes, this is summing booleans

    p_value = n_successes / n_trials

    if verbose:
        print("Monte-Carlo simulation completed. Final p-value: %.6f" % p_value)

    return MonteCarloSimulationResult(p_value, n_trials, n_successes)
