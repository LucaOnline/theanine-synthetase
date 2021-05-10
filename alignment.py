from typing import Tuple


def align_sequences(seq1: str, seq2: str) -> Tuple[str, str]:
    """
    This function aligns the two provided sequences using Needleman-Wunsh
    alignment. It uses a scoring scheme with a gap penalty of -1, a match
    bonus of 1, and a mismatch penalty of 1.
    """
