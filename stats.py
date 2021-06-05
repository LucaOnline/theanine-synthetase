import numpy as np


def variance(data: np.ndarray, sample: bool = True) -> float:
    """
    Calculates the variance within the provided sample or population.
    The `sample` kwarg controls whether the input data should be
    treated as a sample or a population.
    """
    numerator = np.sum(np.square(data) - np.square(data.sum()) / data.size)
    denominator = data.size - 1 if sample else data.size  # if population
    return numerator / denominator
