"""The `aminos` modules contains utilities for interpreting nucleotide sequences in terms of amino acids."""

from typing import Iterator


ALA = "A"
ARG = "R"
ASN = "N"
ASP = "D"
CYS = "C"
GLU = "E"
GLN = "Q"
GLY = "G"
HIS = "H"
ILE = "I"
LEU = "L"
LYS = "K"
MET = "M"
PHE = "F"
PRO = "P"
SER = "S"
THR = "T"
TRP = "W"
TYR = "Y"
VAL = "V"

# Almost all of the mRNA transcripts on NCBI seem to be cDNA,
# so I'm keeping it that way to be faithful to the raw data

STOP_CODONS = ["TAA", "TAG", "TGA"]

AMINO_ACID_TRANSLATIONS = {
    "TTT": PHE,
    "TTC": PHE,
    "TTA": LEU,
    "TTG": LEU,
    "TCT": SER,
    "TCC": SER,
    "TCA": SER,
    "TCG": SER,
    "TAT": TYR,
    "TAC": TYR,
    "TGT": CYS,
    "TGC": CYS,
    "TGG": TRP,
    "CTT": LEU,
    "CTC": LEU,
    "CTA": LEU,
    "CTG": LEU,
    "CCT": PRO,
    "CCC": PRO,
    "CCA": PRO,
    "CCG": PRO,
    "CAT": HIS,
    "CAC": HIS,
    "CAA": GLN,
    "CAG": GLN,
    "CGT": ARG,
    "CGC": ARG,
    "CGA": ARG,
    "CGG": ARG,
    "ATT": ILE,
    "ATC": ILE,
    "ATA": ILE,
    "ATG": MET,
    "ACT": THR,
    "ACC": THR,
    "ACA": THR,
    "ACG": THR,
    "AAT": ASN,
    "AAC": ASN,
    "AAA": LYS,
    "AAG": LYS,
    "AGT": SER,
    "AGC": SER,
    "AGA": ARG,
    "AGG": ARG,
    "GTT": VAL,
    "GTC": VAL,
    "GTA": VAL,
    "GTG": VAL,
    "GCT": ALA,
    "GCC": ALA,
    "GCA": ALA,
    "GCG": ALA,
    "GAT": ASP,
    "GAC": ASP,
    "GAA": GLU,
    "GAG": GLU,
    "GGT": GLY,
    "GGC": GLY,
    "GGA": GLY,
    "GGG": GLY,
}


def changes_amino(codon1: str, codon2: str) -> bool:
    """
    Returns True if codon2's differences from codon1 would change the
    amino acid that codon1 encodes.
    """
    return AMINO_ACID_TRANSLATIONS[codon1] != AMINO_ACID_TRANSLATIONS[codon2]


def codons(sequence: str) -> Iterator[str]:
    """
    Iterates over each codon of the provided sequence. Leading nucleotides
    are discarded until a start codon is reached, and nucleotides after the
    first stop codon are ignored.
    """
    start_index = sequence.find("ATG")
    if start_index == -1:
        raise ValueError("no start codon found")

    sequence = sequence[start_index:]
    codons = [sequence[i * 3 : i * 3 + 3] for i in range(0, len(sequence) // 3)]

    stopped = False
    for codon in codons:
        if codon not in STOP_CODONS:
            yield codon
        else:
            stopped = True
            break

    if not stopped:
        raise ValueError("sequence had no stop codon")
