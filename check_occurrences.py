from typing import Iterator, Tuple

from data_index import CAMELLIA_GENOME, TEA015198
from parse_fasta import parse_fasta


def get_data(filename: str) -> str:
    """Gets the relative file path for the provided file."""
    return "./data/" + filename


def check_occurrences(needle: str, haystacks: Iterator[str]) -> int:
    """
    Returns the number of occurrences of the provided needle
    in the iterator of haystacks.
    """
    return sum(map(lambda haystack: haystack.count(needle), haystacks))


if __name__ == "__main__":
    genome = parse_fasta(get_data(CAMELLIA_GENOME))
    cstsi_mrna = list(parse_fasta(get_data(TEA015198)))[0][1]
    print(
        check_occurrences(
            # Pass in only the actual dataset sequences, not sequence names
            cstsi_mrna,
            map(lambda seq: seq[1], genome),
        )
    )
