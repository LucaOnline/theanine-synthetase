from typing import Callable, TypeVar

T = TypeVar("T")


class MonteCarloSimulationResult:
    """
    MonteCarloSimulationResult represents the result of a Monte-Carlo simulation.
    """

    def __init__(self, p_value: float):
        self.p_value = p_value

    def get_p_value(self):
        """Returns the estimated p-value of the simmulation."""
        return self.p_value


def monte_carlo(
    simulation_fn: Callable[[], T],
    effect_size_fn: Callable[[T], float],
    observed_effect_size: float,
    n_trials: int,
) -> MonteCarloSimulationResult:
    """
    Runs a Monte-Carlo simulation. Returns an object representing the result of the simulation,
    which will contain the final p-value.

    simulation_fn: A function that runs a single trial of the simulation.
    effect_size_fn: A function that calculates an effect size from the return value of the simulation function.
    observed_effect_size: The observed effect size that is being tested.
    n_trials: The number of trials to run. More trials will more precisely estimate the p-value.
    """
    successes = sum(
        [
            1 if effect_size_fn(simulation_fn(trial)) >= observed_effect_size else 0
            for trial in range(n_trials)
        ]
    )
    p_value = successes / n_trials
    return MonteCarloSimulationResult(p_value)
