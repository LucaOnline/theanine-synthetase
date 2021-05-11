from alignment import align_sequences
from data_index import GLUTAMINE_SYNTHETASE, THEANINE_SYNTHETASE
from parse_fasta import parse_fasta


def get_data(filename: str) -> str:
    """Gets the relative file path for the provided file."""
    return "./data/" + filename


if __name__ == "__main__":
    csgsi_mrna = list(parse_fasta(get_data(GLUTAMINE_SYNTHETASE)))[0][1]
    cstsi_mrna = list(parse_fasta(get_data(THEANINE_SYNTHETASE)))[0][1]
    alignment_result = align_sequences("GCATGCT", "GATTACA", nucleotides=True)
    alignment_result.examine()
