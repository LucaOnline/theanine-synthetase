from typing import Dict

from data_index import LONGJING43, TEA015198
from parse_fasta import parse_fasta


def get_data(filename: str) -> str:
    return "./data/" + filename


def check_occurrences(needle: str, haystack: Dict[str, str]) -> int:
    count = 0
    for value in haystack.values():
        if value == needle:
            count += 1
    return count


if __name__ == "__main__":
    longjing43 = parse_fasta(get_data(LONGJING43))
    cstsi = parse_fasta(get_data(TEA015198))
    print(check_occurrences(list(cstsi.values())[0], longjing43))
