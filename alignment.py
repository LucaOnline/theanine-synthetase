from typing import Tuple

import numpy as np


def score_cell(quad: np.ndarray, top_char: str, left_char: str) -> np.int:
    down_score = quad[0, 1] - 1
    right_score = quad[1, 0] - 1
    diag_score = quad[0, 0] - 1
    if top_char == left_char:
        diag_score += 2
    return max([down_score, right_score, diag_score])


def align_sequences(top_seq: str, left_seq: str) -> Tuple[str, str]:
    """
    This function aligns the two provided sequences using Needleman-Wunsh
    alignment. It uses a scoring scheme with a gap penalty of -1, a match
    bonus of 1, and a mismatch penalty of -1.
    """

    size1 = len(top_seq) + 1
    size2 = len(left_seq) + 1

    # Build search matrix
    search = np.zeros((size1, size2), dtype=np.int)
    search[0] = [i for i in range(0, -size1, -1)]
    search[:, 0] = [i for i in range(0, -size2, -1)]

    # Do calculation
    for x in range(1, size1):
        for y in range(1, size2):
            search[x, y] = score_cell(
                search[x - 1 : x + 1, y - 1 : y + 1], top_seq[x - 1], left_seq[y - 1]
            )
    search = search.T

    # Unwind result
    final_top = ""
    final_left = ""

    # TODO backtracking

    print(search)

    return (final_top, final_left)
