import glob
import json
import os
import subprocess

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from biopandas.pdb import PandasPdb

# Helper functions


# read PDB file
def load_pdb_2_df(input_pdb):
    """Loads a pdb file into a Pandas dataframe"""
    df_pdb = PandasPdb().read_pdb(input_pdb).df["ATOM"]
    return df_pdb


# read ColabFold metrics
def load_json_data(target_pdb):
    target_json = get_json_from_pdb(target_pdb)
    with open(target_json, "r") as inpt:
        data_json = json.load(inpt)
    return data_json


def get_json_from_pdb(pdb_name):
    """
    This function will obtain the name of the json file that contains the ColabFold scores
    Only compatible with ColabFold predictions, not ESMfold
    """
    json_file = (
        pdb_name.replace("_unrelaxed", "")
        .replace(".pdb", ".json")
        .replace("_rank", "_scores_rank")
    )
    return json_file


def fix_esmfold_pdb(df_pdb):
    """Fixes the issue with different chains having different residue numbers in ESMfold"""
    df_pdb_fixed = df_pdb.copy()
    # we will apply the index of the first chain to all the other ones
    ref_chain = df_pdb_fixed["chain_id"].unique()[0]
    ref_chain_resi = df_pdb_fixed[df_pdb_fixed["chain_id"] == ref_chain][
        "residue_number"
    ]

    for chain in df_pdb_fixed["chain_id"].unique()[1:]:
        df_pdb_fixed.loc[df_pdb_fixed["chain_id"] == chain, "residue_number"] = (
            ref_chain_resi.values
        )
    df_pdb_fixed["residue_number"] = df_pdb_fixed["residue_number"].astype(int)

    return df_pdb_fixed


def trim_flanking_seq(df_pdb, N_length_flank):
    df_pdb = df_pdb.copy()

    # remove flanking residues
    end_flank = N_length_flank
    start_flank = N_length_flank + fragment_size + 1
    df_pdb_trimmed = df_pdb.loc[
        (df_pdb["residue_number"] > end_flank)
        & (df_pdb["residue_number"] < start_flank),
        :,
    ]

    # reindex atom number
    df_pdb_trimmed.reset_index(inplace=True, drop=True)
    df_pdb_trimmed.loc[:, "atom_number"] = df_pdb_trimmed.index + 1

    # reindex residue number
    df_pdb_trimmed.loc[:, "residue_number"] = (
        df_pdb_trimmed["residue_number"] - df_pdb_trimmed["residue_number"].min() + 1
    )

    return df_pdb_trimmed


def write_updated_pdb(target_pdb, df_pdb_trimmed, output_pdb):
    # write new df_pdb to a pdb file (a little hacky, but it is what it is)
    # load a pandaspdb object for the data structure
    pandaspdb_obj = PandasPdb().read_pdb(target_pdb)
    # substitute the atom section with the new one
    pandaspdb_obj.df["ATOM"] = df_pdb_trimmed
    # write the pdb object to a pdb file
    pandaspdb_obj.to_pdb(output_pdb, records=["ATOM"])


def extract_sequence_from_pandaspdb(df_pdb):
    """Extracts the protein sequence from a PandasPdb"""
    iter_3letter_codes = df_pdb.drop_duplicates(subset=["residue_number"])[
        "residue_name"
    ].values
    prot_sequence = "".join(
        [
            protein_letters_3to1[three_code.capitalize()]
            for three_code in iter_3letter_codes
        ]
    )
    return prot_sequence


# parallel/antiparallel feature


def compute_dist(mat_A, mat_B):
    """Computes the distance between every point in matrix A and matrix B
    The rows are points, the columns coordinates"""
    vec_distances = np.sqrt(np.sum(np.power(mat_A - mat_B, 2), axis=1))
    return vec_distances


def detect_if_pdb_parallel(df_pdb):
    # load CA coordinates
    df_pdb = df_pdb[df_pdb.atom_name == "CA"]

    # get list of chains
    list_chains = df_pdb["chain_id"].unique().tolist()

    # store the coordinates of the CA of every chain
    list_chain_coords = []
    for chain in list_chains:
        coord_chain = df_pdb[df_pdb["chain_id"] == chain][
            ["x_coord", "y_coord", "z_coord"]
        ].values
        list_chain_coords.append(coord_chain)

    # compare distance at each layer between reference chain and any other chain
    ref_coords = list_chain_coords[0]
    result = "PARALLEL"
    for target_coords in list_chain_coords[1:]:
        default_dists = compute_dist(ref_coords, target_coords)
        inverted_dists = compute_dist(ref_coords, target_coords[::-1])
        if np.mean(default_dists) > np.mean(inverted_dists):
            result = "ANTIPARALLEL"

    if result == "PARALLEL":
        parallel = 1
    else:
        parallel = 0

    return parallel


def run_dssp(pdb_file):
    # DSSP4 requires a CRYST1 record, and so it shall have it
    with open(pdb_file, "r") as inpt:
        contents = inpt.readlines()
    if "CRYST1" in contents[0]:
        pass
    else:
        # add the required field at the top of the file
        fix_str = (
            "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"
        )
        contents.insert(0, fix_str)
        with open(pdb_file, "w") as output:
            output.writelines(contents)

    # run DSSP
    dssp_command = f"mkdssp -i {pdb_file} -o temp.dssp"
    subprocess.call(dssp_command, shell=True)
    return 0


def run_socket(pdb_file, threshold=7.4):
    run_dssp(pdb_file)
    # Fix the output to be compatible with Socket
    bad_string = "            CHAIN"
    with open("temp.dssp", "r") as inpt:
        data = inpt.read()
    with open("temp.dssp", "w") as output:
        fixed_data = data.replace(bad_string, "")
        output.write(fixed_data)

    # run socket
    socket_command = f"socket -f {pdb_file} -s temp.dssp -c {threshold} > temp.socket"
    subprocess.call(socket_command, shell=True)
    os.remove("temp.dssp")


def parse_socket(input_socket="temp.socket"):
    with open(input_socket) as inpt:
        list_knobs = [line.strip().split() for line in inpt if line.startswith("knob ")]
    os.remove(input_socket)
    return list_knobs


def detect_kih(input_pdb, length_fragment):
    run_socket(input_pdb)
    list_knobs = parse_socket()
    list_resi_knobs = [data[6].split(":")[0] for data in list_knobs]
    vec_knobs = np.zeros(shape=length_fragment)
    for i in range(1, len(vec_knobs) + 1):
        vec_knobs[i - 1] = list_resi_knobs.count(str(i))
    return vec_knobs


def compute_DSSP_helix(input_pdb):
    # load the input PDB into BioPython parser and run DSSP
    p = PDBParser()
    model = p.get_structure("name_str", input_pdb)[0]
    dssp = DSSP(model, input_pdb, dssp="mkdssp")

    # DSSP data is accessed by a tuple (chain_id, res_id)
    list_SS_symbols = [dssp[i][2] for i in dssp.keys()]

    # Filter for helicity (H)
    list_SS_H = [1 if i == "H" else 0 for i in list_SS_symbols]

    # reshape into nmerxN_residues array
    mat_SS_H = np.array(list_SS_H).reshape(nmer, -1)

    # take as helical residues for which half or more of the chains are helical
    threshold = int(np.ceil(nmer / 2))
    mat_SS_binary = mat_SS_H.sum(axis=0) >= threshold

    return mat_SS_binary


def incorporate_plddt():
    mat_plddt_construct = df_pdb.groupby("residue_number")["b_factor"].mean().values

    mat_results_plddt = dict_result_matrices["plddt"]
    mat_results_plddt[i_construct, row["from"] : row["to"]] = mat_plddt_construct


def incorporate_para():
    mat_para_construct = detect_if_pdb_parallel(df_pdb)

    mat_results_para = dict_result_matrices["para"]
    mat_results_para[i_construct, row["from"] : row["to"]] = mat_para_construct


def incorporate_kih():
    vec_knobs = detect_kih(target_pdb, fragment_size)
    norm_vec_knobs = (
        vec_knobs / nmer
    )  # normalize the number of knobs per position wrt. number of chains

    mat_results_kih = dict_result_matrices["kih"]
    mat_results_kih[i_construct, row["from"] : row["to"]] = norm_vec_knobs


def incorporate_pae():
    pae = np.mean(data_json["pae"])

    mat_results_pae = dict_result_matrices["pae"]
    mat_results_pae[i_construct, row["from"] : row["to"]] = pae


def incorporate_dssp_H():
    mat_heli_construct = compute_DSSP_helix(target_pdb)

    mat_results_heli = dict_result_matrices["heli"]
    mat_results_heli[i_construct, row["from"] : row["to"]] = mat_heli_construct


# this dictionary maps "feature names" to the functions that incorporate them
dict_feature_functions = {
    "plddt": incorporate_plddt,
    "para": incorporate_para,
    "kih": incorporate_kih,
    "pae": incorporate_pae,
    "heli": incorporate_dssp_H,
}


class Integrator:
    def __init__(
        self,
        path_fragments,
        list_features=["plddt"],
        prediction_format="colabfold",
        data_output="DATA_DEIfold.csv",
        model_subfolder="predictions",
        dict_feature_flattening={},
    ):
        """
        prediction_format: colabfold / esmfold
        list_features = ['plddt', 'para', 'kih', 'pae', 'heli']
        """
        self.list_features = list_features
        self.prediction_format = prediction_format.lower()
        self.data_output = data_output
        self.model_subfolder = model_subfolder
        self.root_directory = os.getcwd()
        self.path_fragments = path_fragments
        self.dict_feature_flattening = dict_feature_flattening

        if self.prediction_format == "esmfold" and "pae" in self.list_features:
            print("WARNING PAE is not supported with ESMfold")
            self.list_features.remove("pae")

    def run_assembly_analysis(self):
        # point to parameters
        list_features = self.list_features
        prediction_format = self.prediction_format
        data_output = self.data_output
        model_subfolder = self.model_subfolder

        # change working directory to path where each specification is stored
        os.chdir(self.root_directory)
        os.chdir(self.path_fragments)

        # make relevant variables global to access them from the incorporate_* functions
        global df_pdb
        global row
        global i_construct
        global dict_result_matrices
        global target_pdb
        global fragment_size
        global nmer
        global data_json

        # read full sequence
        source_full_fasta = glob.glob("source_*.fasta")[0]
        with open(source_full_fasta, "r") as inpt:
            record_full_fasta = SeqIO.read(inpt, format="fasta")
        full_sequence = str(record_full_fasta.seq)
        N_length_sequence = len(full_sequence)

        # report detected full sequence and specifications
        print(f"Commencing analysis of {record_full_fasta.id} ...")
        print(
            f"{record_full_fasta.id} read as a sequence of length {N_length_sequence}"
        )

        list_specification_folders = sorted(
            [folder for folder in glob.glob("*/") if not folder.startswith("__")]
        )
        print(f"Found {len(list_specification_folders)} specification folders:")
        for i in list_specification_folders:
            print(f"* {i}")

        # Check whether the analysis was previously run, and if so, just load the data file
        if os.path.exists(self.data_output):
            df_data_sequence = pd.read_csv(self.data_output, index_col=0)
            comp_feats = set(
                [
                    field.split("_")[-1]
                    for field in df_data_sequence.columns
                    if field.count("_") > 1
                ]
            )
            if comp_feats == set(self.list_features):
                print(f"WARNING: Loading previous results {self.data_output}")
                os.chdir(self.root_directory)
                return list_specification_folders, full_sequence, df_data_sequence
            else:
                os.remove(self.data_output)

        # prepare dataframe to store all the DEIfold data
        df_data_sequence = pd.DataFrame([*full_sequence], columns=["residue_name"])
        df_data_sequence["residue_number"] = np.arange(len(full_sequence))

        dict_result_matrices = {}

        for folder_specification in list_specification_folders:
            # load parameters and constructs
            df_constructs = pd.read_csv(
                f"{folder_specification}/constructs.csv", index_col=0
            )
            with open(f"{folder_specification}/parameters.json", "r") as inpt:
                dict_parameters = json.load(inpt)

            # automatically infer the rest of the relevant folders and files
            folder_predictions = f"{folder_specification}/{model_subfolder}"
            list_pdb_files = glob.glob(f"{folder_predictions}/*.pdb")
            # print(len(list_pdb_files))

            N_constructs = df_constructs.shape[0]
            nmer = dict_parameters["nmer"]
            fragment_size = dict_parameters["L"]

            # initialize matrix of results at -1 (later this is used to mask empty values)
            for feature in list_features:
                dict_result_matrices[feature] = (
                    np.zeros(shape=(N_constructs, N_length_sequence)) - 1
                )

            for i_construct, row in df_constructs.iterrows():
                # this dictionary will be passed to every function that computes features
                # it is a filthy solution, but it works
                dict_func_data = {}
                construct_name = row["construct_name"]
                seq_construct = row["sequence"]

                # load relevant data
                if prediction_format == "colabfold":
                    # if we are using ColabFold, take the top-scoring model for this construct
                    target_pdb = [
                        pdb_file
                        for pdb_file in list_pdb_files
                        if construct_name in pdb_file and "rank_001" in pdb_file
                    ][0]
                    df_pdb = load_pdb_2_df(target_pdb)
                    data_json = load_json_data(target_pdb)

                elif prediction_format == "esmfold":
                    # if we are using ESMfold, there will be only one PDB file
                    target_pdb = [
                        pdb_file
                        for pdb_file in list_pdb_files
                        if construct_name in pdb_file
                    ][0]
                    df_pdb = load_pdb_2_df(target_pdb)
                    df_pdb = fix_esmfold_pdb(df_pdb)
                    write_updated_pdb(target_pdb, df_pdb, "temp.pdb")
                    target_pdb = "temp.pdb"

                else:
                    print("Unsupported format")
                    return 0

                # sanity check (i.e., the sequence is what we expected)
                df_pdb = load_pdb_2_df(target_pdb)
                seq_pdb = extract_sequence_from_pandaspdb(df_pdb)
                try:
                    assert seq_construct == seq_pdb
                except AssertionError:
                    print("The sequence and the structure fragments do not match!")
                    print(
                        "This is likely due to a non-valid residue present in the sequence"
                    )
                    print(seq_construct)
                    print(seq_pdb)
                    continue
                    # raise AssertionError

                # if there are any flanking sequence, remove them, and update target_pdb
                try:
                    N_length_flank = len(dict_parameters["flank"])
                    # print(N_length_flank)
                    if N_length_flank > 0:
                        # remove flanking sequences
                        df_pdb = trim_flanking_seq(df_pdb, N_length_flank)
                        # update pdb file
                        write_updated_pdb(target_pdb, df_pdb, "temp.pdb")
                        target_pdb = "temp.pdb"
                except:
                    pass

                for feature in list_features:
                    feature_function = dict_feature_functions[feature]
                    feature_function()

                if prediction_format == "esmfold":
                    os.remove("temp.pdb")

            specification = folder_specification.replace("/", "")
            for feature in list_features:
                mat_results = dict_result_matrices[feature]
                # mask the non-written positions (no overlap)
                mat_results = np.ma.masked_equal(mat_results, value=-1)
                # and write results to dataframe
                if feature in self.dict_feature_flattening.keys():
                    flattening_function = self.dict_feature_flattening[feature]
                else:
                    flattening_function = np.mean

                df_data_sequence[f"{specification}_{feature}"] = flattening_function(
                    mat_results, axis=0
                )

        # save data to file
        df_data_sequence.to_csv(data_output)

        # return to original working directory
        os.chdir(self.root_directory)

        return list_specification_folders, full_sequence, df_data_sequence
