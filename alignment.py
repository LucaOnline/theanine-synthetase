from typing import Tuple, Literal

import numpy as np


MOVE_DIAGONAL = 0
MOVE_RIGHT = 1
MOVE_DOWN = 2
EditMove = Literal[MOVE_DIAGONAL, MOVE_RIGHT, MOVE_DOWN]


def backtrack(quad: np.ndarray) -> EditMove:
    """Trace one step back through an edit matrix."""
    if quad.shape == (0, 2):
        return MOVE_DOWN
    elif quad.shape == (2, 0):
        return MOVE_RIGHT

    next_pos = (0, 0)
    if quad[0, 1] > quad[next_pos]:
        next_pos = (0, 1)
    if quad[1, 0] > quad[next_pos]:
        next_pos = (1, 0)

    if next_pos == (0, 0):
        return MOVE_DIAGONAL
    elif next_pos == (0, 1):
        return MOVE_RIGHT
    else:
        return MOVE_DOWN


def score_cell(quad: np.ndarray, top_char: str, left_char: str) -> np.int:
    """Calculate the Needleman-Wunsch score for a cell."""
    down_score = quad[0, 1] - 1
    right_score = quad[1, 0] - 1
    diag_score = quad[0, 0] - 1
    if top_char == left_char:
        diag_score += 2
    return max([down_score, right_score, diag_score])


def align_sequences(top_seq: str, left_seq: str) -> Tuple[str, str]:
    """
    This function aligns the two provided sequences using Needleman-Wunsch
    alignment. It uses a scoring scheme with a gap penalty of -1, a match
    bonus of 1, and a mismatch penalty of -1.
    """

    size1 = len(top_seq) + 1
    size2 = len(left_seq) + 1

    # Build search matrix
    search = np.zeros((size1, size2), dtype=np.int)
    search[0] = [i for i in range(0, -size1, -1)]
    search[:, 0] = [i for i in range(0, -size2, -1)]

    # Do scoring
    for x in range(1, size1):
        for y in range(1, size2):
            search[x, y] = score_cell(
                search[x - 1 : x + 1, y - 1 : y + 1], top_seq[x - 1], left_seq[y - 1]
            )
    search = search.T

    # Unwind result
    final_top = ""
    final_left = ""

    bt_x, bt_y = (size1 - 1, size2 - 1)
    while bt_x != 0 or bt_y != 0:
        next_move = backtrack(search[bt_x - 1 : bt_x + 1, bt_y - 1 : bt_y + 1])
        if next_move == MOVE_DIAGONAL:
            final_top = top_seq[bt_x - 1] + final_top
            final_left = left_seq[bt_y - 1] + final_left
            bt_x -= 1
            bt_y -= 1
        elif next_move == MOVE_DOWN:
            final_top = "-" + final_top
            final_left = left_seq[bt_y - 1] + final_left
            bt_y -= 1
        elif next_move == MOVE_RIGHT:
            final_top = top_seq[bt_x - 1] + final_top
            final_left = "-" + final_left
            bt_x -= 1

    return (final_top, final_left)
