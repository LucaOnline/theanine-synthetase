from typing import Iterator, Tuple

from data_index import CAMELLIA_GENOME, TEA015198
from parse_fasta import parse_fasta


def get_data(filename: str) -> str:
    return "./data/" + filename


def check_occurrences(needle: str, haystacks: Iterator[Tuple[str, str]]) -> int:
    return list(map(lambda haystack: 1 if needle in haystack[1] else 0, haystacks))


if __name__ == "__main__":
    genome = parse_fasta(get_data(CAMELLIA_GENOME))
    cstsi_mrna = list(parse_fasta(get_data(TEA015198)))[0][1]
    print(check_occurrences(cstsi_mrna, genome))
