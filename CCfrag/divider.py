import json
import os

import numpy as np
import pandas as pd

from .common import read_fasta


class Divider:
    def __init__(self, L, O, nmer, flank=None, tag=None):
        self.L = L
        self.O = O
        self.nmer = nmer
        self.flank = flank
        self.tag = tag

    def _generate_fragment_idx(self, sequence):
        """
        This function generates pairs of indices that divide the input sequence into the
        specification defined by the user
        """
        seq_len = len(sequence)

        # if the user specified no fragments
        if self.L == -1:
            list_segment_idx = [(0, seq_len)]

        # if the user specified fragment size
        else:
            # compute the indices for each fragment
            list_segment_idx = [
                (i, i + self.L) for i in list(range(0, len(sequence), self.L - self.O))
            ]

            # remove any fragments at the end that are not complete
            list_segment_idx = list(filter(lambda x: x[1] <= seq_len, list_segment_idx))

            # add the last complete fragment (might overlap by more than O with previous fragment)
            last_idx = (seq_len - self.L, seq_len)
            list_segment_idx.append(last_idx)

            # check whether the two last indices are identical
            # (can happen if overlap and window sizes align)
            if list_segment_idx[-2] == list_segment_idx[-1]:
                list_segment_idx = list_segment_idx[:-1]

        return list_segment_idx

    def generate_queries(
        self, input_fasta, output_path, format="single_fasta", overwrite=False
    ):
        """
        This is the top-level operating function of the disassembler
        """
        # for each sequence in the input_fasta file, apply the windowing
        for main_header, full_sequence in read_fasta(input_fasta):
            print(main_header)

            # determine the name of the sequence
            seq_name = main_header.split()[0].replace(">", "")

            # sanitize input parameters (check for errors)
            self._sanitize_parameters(full_sequence)

            # compute list of indices
            list_segment_idx = self._generate_fragment_idx(full_sequence)

            # compute list of subsequences
            list_segments = [full_sequence[i:j] for i, j in list_segment_idx]

            # if provided, add the flanking sequences
            if self.flank is not None and len(self.flank) > 0:
                list_segments = [self.flank + seq + self.flank for seq in list_segments]

            # compute list of queries
            N = len(list_segments)
            list_queries = list(
                zip(
                    range(N),
                    list_segments,
                    list_segment_idx,
                    [self.nmer] * N,
                    [seq_name] * N,
                )
            )

            # export queries to csv format
            # create queries dataframe
            df_queries = pd.DataFrame(
                list_queries, columns=["i", "sequence", "idx", "nmer", "seq_name"]
            )
            df_queries["from"] = df_queries["idx"].apply(lambda x: x[0])
            df_queries["to"] = df_queries["idx"].apply(lambda x: x[1])
            df_queries["construct_name"] = (
                df_queries["seq_name"]
                + "_"
                + df_queries["nmer"].astype(str)
                + "_"
                + df_queries["from"].astype(str)
                + "-"
                + df_queries["to"].astype(str)
            )

            # create sequence-specific folder (to store all sequence-derived constructs)
            sequence_output_path = f"{output_path}/{seq_name}"
            if self.tag is not None:
                sequence_output_path = f"{output_path}/{seq_name}_{self.tag}"
            os.makedirs(sequence_output_path, exist_ok=True)

            # create specification-specific folder
            specification_name = f"{self.nmer}_{self.L}_{self.O}"
            specification_path = f"{sequence_output_path}/{specification_name}"
            os.makedirs(specification_path, exist_ok=overwrite)
            df_queries.to_csv(f"{specification_path}/constructs.csv")

            # create the actual queries
            os.makedirs(f"{specification_path}/queries", exist_ok=overwrite)
            for i, row in df_queries.iterrows():
                name = f"{row.construct_name}"
                header = f">{row.construct_name}"
                sequence = f"{row.sequence}"

                if format == "single_fasta":
                    with open(
                        f"{specification_path}/queries/{name}.fasta", "w"
                    ) as output:
                        output.write(header + "\n")
                        output.write(":".join([sequence] * self.nmer))
                elif format == "multi_fasta":
                    with open(
                        f"{specification_path}/queries/{name}.fasta", "w"
                    ) as output:
                        for j in range(self.nmer):
                            output.write(header + "\n")
                            output.write(sequence + "\n")
                else:
                    print("ERROR! Please use format 'single_fasta' or 'multi_fasta'")
                    return 0

            # export parameters to a csv format
            dict_params = {
                "nmer": self.nmer,
                "L": self.L,
                "O": self.O,
                "flank": self.flank,
            }
            with open(f"{specification_path}/parameters.json", "w") as output:
                json.dump(dict_params, output)

            # copy source sequence to the output path (to be used in the analysis part)
            output_source_file = f"{sequence_output_path}/source_{seq_name}.fasta"
            with open(output_source_file, "w") as output:
                output.write(f"{main_header}\n")
                output.write(f"{full_sequence}\n")

    def _sanitize_parameters(self, sequence):
        """
        This function checks that the input parameters are feasible and that will not be conducive to error
        """
        # The sequence must consist solely of standard amino acids in their one-letter format.
        standard_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")  # Standard 20 amino acids
        if all(residue in standard_amino_acids for residue in sequence.upper()):
            pass
        else:
            errmsg = "Non-standard residues detected in input sequence!"
            raise Exception(errmsg)

        # L must be either -1 (no fragment modeling), or a minimum of 10
        if self.L == -1 or self.L >= 10:
            pass
        else:
            errmsg = "L (window_size) must either be -1 (no fragment modeling), or a minimum of 10"
            raise Exception(errmsg)

        if self.L == -1:
            return 0

        # Window size must be smaller than sequence length
        if self.L < len(sequence):
            pass
        else:
            errmsg = "L (window_size) must be smaller than sequence_length!"
            raise Exception(errmsg)

        # O (overlap_size) must be between 0 (no overlap) and window_length - 1 (maximum possible overlap)
        if self.O >= 0 and self.O < self.L:
            pass
        else:
            errmsg = "O (overlap_length) must be smaller than window_length!"
            raise Exception(errmsg)

    def debug_check(self):
        """
        This function will check that
        * The fragment generation function:
            + Generates all fragments of equal size
            + Generates fragments that cover the entirety of the sequence
        """
        # generate a test sequence
        seqAA = "".join(list("ACDEFGHIKL"))
        sequence = seqAA * 4 + seqAA[::-1]

        print(f"Testing fragment generation with L={self.L} and O={self.O}")
        self._sanitize_parameters(sequence)

        # visually check that size of all fragments are identical
        list_segment_idx = self._generate_fragment_idx(sequence)
        print(list_segment_idx)

        # check that the fragments cover 100% of the sequence
        vec_test = np.zeros(len(sequence))

        for i, j in list_segment_idx:
            segment = sequence[i:j]
            print(" " * i + segment)

            # check that each fragment is of the correct size
            assert len(segment) == self.L
            vec_test[i:j] = 1

        assert vec_test.sum() == len(sequence)

        print(f"No errors found with L={self.L} and O={self.O}")
        print()
