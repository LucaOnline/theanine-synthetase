from itertools import product


NUCLEOTIDE_ORDER = ["T", "C", "A", "G"]
AMINO_ACIDS = [
    "F",
    "L",
    "S",
    "Y",
    "!", # Sentinel for STOP
    "C",
    "!", # Sentinel for STOP
    "W",
    "L",
    "P",
    "H",
    "Q",
    "R",
    "I",
    "M",
    "T",
    "N",
    "K",
    "S",
    "R",
    "V",
    "A",
    "D",
    "E",
    "G",
]

def get_amino_acid_table():
    """
    Returns a table keyed with nucleotide sequences and their translations in
    terms of amino acids.
    """
    aminos = ["".join(tri) for tri in product("".join(NUCLEOTIDE_ORDER), repeat=3)]
    aminos = sorted(
        aminos,
        key=lambda amino_acid: [NUCLEOTIDE_ORDER.index(b) for b in amino_acid]
    )
    amino_groups = [aminos[i:i+4] for i in range(0, len(aminos), 4)]

    temp_amino_groups = []
    for group in amino_groups:
        if group[0][1] == "A" or group[0][0:2] == "TT" or group[0][0:2] == "AG":
            temp_amino_groups.append(group[0:2])
            temp_amino_groups.append(group[2:4])
        elif group[3] == "ATG":
            temp_amino_groups.append(group[0:3])
            temp_amino_groups.append([group[3]])
        elif group[2] == "TGA":
            temp_amino_groups.append(group[0:2])
            temp_amino_groups.append([group[2]])
            temp_amino_groups.append([group[3]])
        else:
            temp_amino_groups.append(group)
    amino_groups = temp_amino_groups

    amino_table = {}
    for group, amino in zip(amino_groups, AMINO_ACIDS):
        for seq in group:
            amino_table[seq] = amino

    return amino_table

AMINO_TABLE = get_amino_acid_table()


def to_protein(nucleotides: str) -> str:
    """
    Converts a nucleotide sequence into a single-letter code amino acid
    sequence. Triplet nucleotides are taken from the provided input until
    this is no longer possible.
    """
    return [AMINO_TABLE[nucleotides[i:i+3]] for i in range(0, (len(nucleotides) // 3) * 3, 3)]