from typing import Iterator, Tuple


def parse_fasta(filename: str) -> Iterator[Tuple[str, str]]:
    with open(filename, "r") as f:
        data = f.read()

    if data[0] != ">":
        raise ValueError("input file is not a FASTA file")

    currentSeq = []
    currentKey = ""
    for line in data.splitlines():
        if line.startswith(">"):
            if currentKey != "":
                yield currentKey, "".join(currentSeq)
                currentSeq = []
            currentKey = line[1:].split()[0]  # Only key on the protein name
        else:
            currentSeq.append(line.strip().upper())

    # Return last pair, if it exists
    if currentKey != "":
        yield currentKey, "".join(currentSeq)
