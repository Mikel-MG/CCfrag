# fasta file parsing function (multi-line capable)
def read_fasta(filename):
    """
    This function parses a FASTA file and returns a generator that iterates through each header and sequence
    """
    with open(filename, 'r') as inpt:
        header = None
        sequence = []

        for line in inpt:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence).upper()
                header = line
                sequence = []
            else:
                sequence.append(line)

        if header is not None:
            yield header, ''.join(sequence).upper()
