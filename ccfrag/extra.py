import glob
import json
import pickle
import shutil
from pathlib import Path


def extract_best_model(ranking_json: Path):
    """
    Extracts the name of the best-ranking model from a json file
    """
    with open(ranking_json, "r") as inpt:
        data_ranking = json.load(inpt)
    best_model = data_ranking["order"][0]
    return best_model


def format_scores(input_pkl: Path, output_json: Path):
    """
    Converts the prediction confidence metrics from a .pkl file to a .json file
    """
    # load the pkl file data
    with open(input_pkl, "rb") as inpt:
        data = pickle.load(inpt)

    # reformat to match the colabfold style
    data["pae"] = data["predicted_aligned_error"].tolist()
    data.pop("predicted_aligned_error")
    data["plddt"] = data["plddt"].tolist()
    data["ptm"] = data["ptm"].tolist()
    data["iptm"] = data["iptm"].tolist()

    # save score file as a json
    with open(output_json, "w") as output:
        json.dump(data, output)


class AF2_to_ColabFold_formatter:
    def __init__(
        self,
        path_input_AF2_results: Path,
        path_output_ColabFold_results: Path,
        prediction_folder_name: str = "predictions",
    ):
        """
        relpath_AF2_results/
         * 2_20_10
         * 2_30_15
         source_SEQ.fasta
         ...
        """
        self.path_input_AF2_results = path_input_AF2_results
        self.path_output_ColabFold_results = path_output_ColabFold_results
        self.prediction_folder_name = prediction_folder_name

    def reformat_results(self):
        # make folder for the reformatted results
        # note: it will remove old reformatted results
        shutil.rmtree(self.path_output_ColabFold_results, ignore_errors=True)
        self.path_output_ColabFold_results.mkdir()

        # copy source fasta file
        source_fasta = list(self.path_input_AF2_results.glob("source*.fasta"))[0]
        shutil.copy(source_fasta, self.path_output_ColabFold_results)

        # list_specification_folders = glob.glob(f"{self.relpath_AF2_results}/*/")
        list_specification_folders = list(self.path_input_AF2_results.glob("*/"))

        for specification_path in list_specification_folders:
            # identify the specification
            specification_name = specification_path.stem

            # create folder for formatted output
            path_specification_output = (
                self.path_output_ColabFold_results / specification_name
            )
            path_specification_output.mkdir()

            path_prediction_folder = (
                path_specification_output / self.prediction_folder_name
            )
            path_prediction_folder.mkdir()

            # copy parameter and construct files
            parameters_file = f"{specification_path}/parameters.json"
            constructs_file = f"{specification_path}/constructs.csv"

            shutil.copy(parameters_file, path_specification_output)
            shutil.copy(constructs_file, path_specification_output)

            # list folders with AlphaFold predictions
            list_prediction_folders = [
                i.parent
                for i in self.path_input_AF2_results.rglob("ranking_debug.json")
            ]
            print(f"Converting {len(list_prediction_folders)} predictions...")

            # go through all the AF2 folders and export them into a common ColabFold folder
            for path_AF2_files in list_prediction_folders:
                # extract best model
                ranking_json = path_AF2_files / "ranking_debug.json"
                best_model = extract_best_model(ranking_json)

                # find best model and corresponding pkl file
                best_pdb = glob.glob(f"{path_AF2_files}/*{best_model}*.pdb")[0]
                best_pkl = glob.glob(f"{path_AF2_files}/*{best_model}*.pkl")[0]

                # generate name of new files
                folder = path_AF2_files.stem
                name_new_pdb = f"{folder}_unrelaxed_rank_001.pdb"
                name_new_json = f"{folder}_scores_rank_001.json"

                # copy pdb file of best-ranking model
                path_new_pdb = f"{path_specification_output}/{self.prediction_folder_name}/{name_new_pdb}"
                shutil.copy(best_pdb, path_new_pdb)

                # produce new score file (in json format)
                path_new_json = path_prediction_folder / name_new_json
                format_scores(Path(best_pkl), path_new_json)
