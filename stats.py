from typing import List


def variance(data: List[int], sample: bool = True) -> float:
    """
    Calculates the variance within the provided sample or population.
    The `sample` kwarg controls whether the input data should be
    treated as a sample or a population.
    """
    numerator = sum(map(lambda x: x * x, data)) - ((sum(data) ** 2) / len(data))
    denominator = len(data) - 1 if sample else len(data)  # if population
    return numerator / denominator
