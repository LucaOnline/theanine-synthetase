from typing import Dict


def parse_fasta(filename: str) -> Dict[str, str]:
    with open(filename, "r") as f:
        data = f.read()

    sequences = {}
    currentKey = ""
    for line in data.splitlines():
        if line.startswith(">"):
            currentKey = line[1:].split("\t")[0]  # Only key on the protein name
        else:
            if sequences.get(currentKey) is None:
                sequences[currentKey] = ""
            sequences[currentKey] += line.strip().upper()
    return sequences
