"""
The `dnds_prep` module exposes methods for preparing an alignment for dN/dS (Ka/Ks) analysis.
"""

from typing import Tuple

from alignment import AlignmentResult


def trim_indels(alignment: AlignmentResult) -> Tuple[str, str]:
    """
    Trims indels from the provided alignment and returns a
    tuple (alignment_1, alignment_2).
    """
    alignment_1 = alignment.get_alignment_1()
    alignment_2 = alignment.get_alignment_2()
    return alignment_1.replace("-", ""), alignment_2.replace("-", "")
