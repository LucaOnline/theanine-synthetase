"""
The `dnds_prep` module exposes methods for preparing an alignment for dN/dS (Ka/Ks) analysis.
"""

from typing import Tuple

from alignment import AlignmentResult


def trim_for_dnds(alignment: AlignmentResult) -> Tuple[str, str]:
    """
    Trims indels from the provided alignment, along with trailing nucleotides,
    and returns a tuple (alignment_1, alignment_2).
    """

    alignment_1 = alignment.get_alignment_1()
    alignment_2 = alignment.get_alignment_2()

    indels_1 = [i for i, c in enumerate(alignment_1) if c == "-"]
    indels_2 = [i for i, c in enumerate(alignment_2) if c == "-"]
    indels_all = set(indels_1 + indels_2)

    trimmed_1 = "".join([c for i, c in enumerate(alignment_1) if not i in indels_all])
    trimmed_2 = "".join([c for i, c in enumerate(alignment_2) if not i in indels_all])

    trimmed_1 = "".join(
        [trimmed_1[i * 3 : i * 3 + 3] for i in range(0, len(trimmed_1) // 3)]
    )
    trimmed_2 = "".join(
        [trimmed_2[i * 3 : i * 3 + 3] for i in range(0, len(trimmed_2) // 3)]
    )

    return trimmed_1, trimmed_2
