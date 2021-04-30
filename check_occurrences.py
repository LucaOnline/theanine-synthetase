from typing import Iterator, Tuple
from multiprocessing import Pool

from data_index import CAMELLIA_GENOME, TEA015198
from parse_fasta import parse_fasta


NUCLEOTIDE_INVERSIONS = {
    "G": "C",
    "C": "G",
    "A": "T",
    "T": "A",
}


def get_data(filename: str) -> str:
    return "./data/" + filename


def invert_nucleotides(nucleotides: str) -> str:
    return "".join([NUCLEOTIDE_INVERSIONS[n] for n in nucleotides])


def check_occurrences(needle: str, haystacks: Iterator[Tuple[str, str]]) -> int:
    with Pool() as pool:
        return len(pool.map(lambda haystack: haystack[1].count(needle), haystacks))


if __name__ == "__main__":
    genome = parse_fasta(get_data(CAMELLIA_GENOME))
    cstsi_mrna = list(parse_fasta(get_data(TEA015198)))[0][1]
    print(check_occurrences(invert_nucleotides(cstsi_mrna), genome))
