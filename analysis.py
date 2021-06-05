import json

from alignment import align_sequences
from data_index import CSTSI, CSGSI
from dir_utils import get_data, make_output_dir, get_output
from options import CLUSTER_COUNT
from parse_fasta import parse_fasta


if __name__ == "__main__":
    # Read CsTSI sequence
    cstsi_mrna = list(parse_fasta(get_data(CSTSI)))[0][1]

    # Analyze CsGSI sequence
    print("Analyzing %s..." % CSGSI)

    csgsi_seq = list(parse_fasta(get_data(CSGSI)))[0][1]
    alignment_result = align_sequences(csgsi_seq, cstsi_mrna, nucleotides=True)

    largest_mismatch_pos, largest_mismatch = alignment_result.largest_mismatch()
    percent_similarity = 1 - (
        alignment_result.hamming_distance() / alignment_result.get_alignment_length()
    )

    clustered_mismatches = alignment_result.clustered_mismatches(
        cluster_count=CLUSTER_COUNT
    )
    clustered_mismatch_variance = alignment_result.clustered_mismatch_variance(
        cluster_count=CLUSTER_COUNT
    )

    # Output Supplementary Data 4
    make_output_dir()
    with open(get_output(CSGSI + ".aln.txt"), "w+") as f:
        f.write(alignment_result.format_result(line_length=100))

    with open(get_output(CSGSI + ".meta.txt"), "w+") as f:
        f.write(
            "Formatted metadata -- not for programmatic use.\n\n"
            + "Information for alignment with CsTSI:\n\n"
            + f"Percent similarity: {percent_similarity}\n"
            + f"Largest mismatch location: {largest_mismatch_pos}\n"
            + f"Largest mismatch size: {largest_mismatch}bp\n"
            + f"Variance between clusters ({CLUSTER_COUNT} clusters): {clustered_mismatch_variance}\n"
            + f"Clustered mismatches: {clustered_mismatches}\n"
        )

    with open(get_output(CSGSI + ".meta.json"), "w+") as f:
        json.dump(
            {
                "percent_similarity": percent_similarity,
                "largest_mismatch_pos": largest_mismatch_pos,
                "largest_mismatch": largest_mismatch,
                "clustered_mismatch_variance": clustered_mismatch_variance,
                "clustered_mismatches": clustered_mismatches,
            },
            f,
        )
