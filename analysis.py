import os

from alignment import align_sequences
from data_index import DATA_URLS, CSTSI
from parse_fasta import parse_fasta


def get_data(filename: str) -> str:
    """Gets the relative file path for the provided file."""
    return "./data/" + filename


def make_output_dir():
    """Creates the output file directory, if it doesn't exist."""
    if not os.path.exists("./output"):
        os.mkdir("./output")


def get_output(filename: str) -> str:
    """Gets the relative file path for the provided output file."""
    return "./output/" + filename


if __name__ == "__main__":
    cstsi_mrna = list(parse_fasta(get_data(CSTSI)))[0][1]

    # Analyze all GS sequences downloaded
    test_sequences = [seq for seq in DATA_URLS.keys() if seq != CSTSI]
    for sequence in test_sequences:
        print("Analyzing %s..." % sequence)

        test_seq = list(parse_fasta(get_data(sequence)))[0][1]
        alignment_result = align_sequences(test_seq, cstsi_mrna, nucleotides=True)

        print(alignment_result.clustered_mismatches(15))

        largest_mismatch_pos, largest_mismatch = alignment_result.largest_mismatch()
        percent_similarity = 1 - (
            alignment_result.hamming_distance()
            / alignment_result.get_alignment_length()
        )

        # Output Supplementary Data 4
        make_output_dir()
        with open(get_output(sequence + ".aln.txt"), "w+") as f:
            f.write(alignment_result.format_result(line_length=100))
        with open(get_output(sequence + ".meta.txt"), "w+") as f:
            f.write(
                f"Percent similarity: {percent_similarity}\nLargest mismatch location: {largest_mismatch_pos}\nLargest mismatch size: {largest_mismatch}bp\n"
            )
