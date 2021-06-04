from typing import Tuple, Literal, List
from math import ceil

import numpy as np

from stats import variance


MOVE_DIAGONAL = 0
MOVE_RIGHT = 1
MOVE_DOWN = 2
EditMove = Literal[MOVE_DIAGONAL, MOVE_RIGHT, MOVE_DOWN]


CHEMICAL_CLASS = {
    "A": "Purine",
    "G": "Purine",
    "T": "Pyrimidine",
    "C": "Pyrimidine",
}


class AlignmentResult:
    """
    AlignmentResult represents the result of performing an alignment on
    two sequences.
    """

    def __init__(self, alignment_1: str, alignment_2: str):
        if len(alignment_1) != len(alignment_2):
            raise ValueError("input strings have differing lengths")
        self.alignment_1 = alignment_1
        self.alignment_2 = alignment_2

    def get_alignment_length(self) -> int:
        """Returns the length of the alignment."""
        return len(self.alignment_1)

    def get_alignment_1(self) -> str:
        """Returns the first alignment string."""
        return self.alignment_1

    def get_alignment_2(self) -> str:
        """Returns the second alignment string."""
        return self.alignment_2

    def get_match_string(self) -> str:
        """Returns the match string for the alignment."""
        return "".join(
            [
                "|" if self.alignment_1[i] == self.alignment_2[i] else " "
                for i in range(len(self.alignment_1))
            ]
        )

    def clustered_mismatches(self, cluster_count: int = 6) -> List[int]:
        """
        Breaks the alignment into `cluster_count` clusters and
        returns the number of mismatches in each cluster. If the
        alignment cannot be equally divided into the number of
        clusters, this leaves the last cluster with the remainder
        of the mismatches.
        """
        if cluster_count < 1:
            raise ValueError("cluster count must be greater than or equal to 1")

        match_string = self.get_match_string()

        cluster_size = ceil(len(match_string) / cluster_count)

        return [
            match_string[i * cluster_size : i * cluster_size + cluster_size].count(" ")
            for i in range(0, len(match_string) // cluster_size)
        ]

    def clustered_mismatch_variance(self, cluster_count: int = 6) -> float:
        """
        Returns the variance between the mismatch clusters. The
        raw cluster mismatches can be retrieved with the
        `clustered_mismatches` method. `cluster_count` controls
        the number of clusters used.
        """
        return variance(
            self.clustered_mismatches(cluster_count=cluster_count), sample=False
        )

    def matches(self) -> int:
        """Returns the number of matching elements for the alignment."""
        return self.get_match_string().count("|")

    def hamming_distance(self) -> int:
        """Returns the Hamming distance of the alignment."""
        return len(self.alignment_1) - self.matches()

    def largest_mismatch(self) -> Tuple[int, int]:
        """Returns the position and size of the largest mismatch in the alignment."""
        matches = self.get_match_string()
        found_mismatch = False
        largest_mismatch = 0
        largest_mismatch_pos = 0
        current_mismatch = 0
        for i, c in enumerate(matches):
            if c == " ":
                found_mismatch = True
                current_mismatch += 1
                if current_mismatch > largest_mismatch:
                    largest_mismatch = current_mismatch
                    largest_mismatch_pos = i - largest_mismatch + 1
            else:
                current_mismatch = 0
        if found_mismatch:
            return (largest_mismatch_pos, largest_mismatch)
        return (-1, 0)

    def format_result(self, line_length: int = 80):
        """
        Formats the found alignment with pipes between
        matching elements. The optional `line_length` parameter
        allows for adjusting the number of elements on each set of
        lines.
        """
        matches = self.get_match_string()

        # Chunk lines
        alignment_1_lines = [
            self.alignment_1[i : i + line_length]
            for i in range(0, len(self.alignment_1), line_length)
        ]
        alignment_2_lines = [
            self.alignment_2[i : i + line_length]
            for i in range(0, len(self.alignment_2), line_length)
        ]
        match_lines = [
            matches[i : i + line_length] for i in range(0, len(matches), line_length)
        ]

        # Output line chunks in order
        return "\n".join(
            [
                "\n".join(
                    [alignment_1_lines[i], match_lines[i], alignment_2_lines[i], ""]
                )
                for i in range(len(match_lines))
            ]
        )

    def examine(self, line_length: int = 80):
        """
        Formats and prints the found alignment with pipes between
        matching elements. The optional `line_length` parameter
        allows for adjusting the number of elements on each set of
        lines.
        """
        print(self.format_result(line_length=line_length))


def backtrack(quad: np.ndarray) -> EditMove:
    """Trace one step back through an edit matrix."""
    if quad.shape == (0, 2):
        return MOVE_DOWN
    elif quad.shape == (2, 0):
        return MOVE_RIGHT

    # numpy's argmax doesn't allow for prioritizing non-indels
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


def score_cell(
    quad: np.ndarray, top_char: str, left_char: str, nucleotides: bool
) -> np.int:
    """Calculate the Needleman-Wunsch score for a cell."""
    down_score = quad[0, 1] - 1
    right_score = quad[1, 0] - 1

    # Penalize transversions more heavily
    if nucleotides and CHEMICAL_CLASS[top_char] != CHEMICAL_CLASS[left_char]:
        down_score -= 1
        right_score -= 1

    diag_score = quad[0, 0] - 1
    if top_char == left_char:
        diag_score += 2
    return max([down_score, right_score, diag_score])


def align_sequences(
    top_seq: str, left_seq: str, nucleotides: bool = True
) -> AlignmentResult:
    """
    This function aligns the two provided sequences using Needleman-Wunsch
    alignment. It uses a scoring scheme with a gap penalty of -1, a match
    bonus of 1, and a mismatch penalty of -1. If the two sequences are
    `nucleotides`, then an additional -1 penalty is applied to transversions.
    """

    size1 = len(top_seq) + 1
    size2 = len(left_seq) + 1

    # Build search matrix
    search = np.zeros((size2, size1), dtype=np.int)
    search[0] = [i for i in range(0, -size1, -1)]
    search[:, 0] = [i for i in range(0, -size2, -1)]

    # Do scoring
    for x in range(1, size2):
        for y in range(1, size1):
            search[x, y] = score_cell(
                search[x - 1 : x + 1, y - 1 : y + 1],
                top_seq[y - 1],
                left_seq[x - 1],
                nucleotides,
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

    return AlignmentResult(final_top, final_left)
