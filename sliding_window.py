"""
The `sliding_window` module provides a function that generates a sliding-window iterable
over a sequence.
"""

from itertools import islice
from typing import Iterator


def sliding_window(sequence: str, n: int) -> Iterator[str]:
    """
    Returns a sliding window (of width n) over data from the iterable,
    e.g. s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...

    Ripped from the Python documentation:
    https://docs.python.org/release/2.3.5/lib/itertools-example.html
    """
    it = iter(sequence)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result
