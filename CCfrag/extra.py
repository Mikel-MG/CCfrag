import os
import glob
import shutil

import json

def extract_best_model(ranking_json):
    with open(ranking_json, 'r') as inpt:
        data_ranking = json.load(inpt)
    best_model = data_ranking['order'][0]
    return best_model

import pickle

def format_scores(input_pkl, output_json):
    # load the pkl file data
    with open(input_pkl, 'rb') as inpt:
        data = pickle.load(inpt)
        
    # reformat to match the colabfold style
    data['pae'] = data['predicted_aligned_error'].tolist()
    data['plddt'] = data['plddt'].tolist()
    data['ptm'] = data['ptm'].tolist()
    data['iptm'] = data['iptm'].tolist()
    data.pop('predicted_aligned_error')
    
    # save score file as a json
    with open(output_json, 'w') as output:
        json.dump(data, output)


class AF2_to_ColabFold_formatter():
    def __init__(self, relpath_AF2_results, relpath_ColabFold_results,
                 prediction_folder_name='predictions'):
        """
        relpath_AF2_results/
         * 2_20_10
         * 2_30_15
         source_SEQ.fasta
         ...
        """
        self.relpath_AF2_results = relpath_AF2_results
        self.relpath_ColabFold_results = relpath_ColabFold_results
        self.prediction_folder_name = prediction_folder_name
        
    def reformat_results(self):
        # make folder for the reformatted results
        # note: it will remove old reformatted results
        shutil.rmtree(self.relpath_ColabFold_results, ignore_errors=True)
        os.makedirs(self.relpath_ColabFold_results)

        # copy source fasta file
        source_fasta = glob.glob(f'{self.relpath_AF2_results}/source*.fasta')[0]
        shutil.copy(source_fasta, self.relpath_ColabFold_results)
        
        # list specifications
        list_specification_folders = glob.glob(f'{self.relpath_AF2_results}/*/')
        
        for specification_path in list_specification_folders:
            # identify the specification
            _, specification_name = os.path.split(os.path.dirname(specification_path))

            # create folder for formatted output
            specification_output = f'{self.relpath_ColabFold_results}/{specification_name}'
            os.makedirs(specification_output)
            os.makedirs(f'{specification_output}/{self.prediction_folder_name}')

            # copy parameter and construct files
            parameters_file = f'{specification_path}/parameters.json'
            constructs_file = f'{specification_path}/constructs.csv'

            shutil.copy(parameters_file, specification_output)
            shutil.copy(constructs_file, specification_output)

            # list folders with AlphaFold predictions
            list_prediction_folders = [i.replace('/ranking_debug.json','').replace('//','/') for i in
                                       glob.glob(f'{self.relpath_AF2_results}/{specification_name}/*/*/ranking_debug.json')]
            list_prediction_folders = sorted(list_prediction_folders, key=lambda x: int(x.split('-')[-1]))

            # Go through all the AF2 folders and export them into a common ColabFold folder
            for path in list_prediction_folders:
                # extract best model
                ranking_json = f'{path}/ranking_debug.json'
                best_model = extract_best_model(ranking_json)

                # find best model and corresponding pkl file
                best_pdb = glob.glob(f'{path}/*{best_model}*.pdb')[0] 
                best_pkl = glob.glob(f'{path}/*{best_model}*.pkl')[0]

                # generate name of new files
                folder = os.path.basename(path)
                name_new_pdb = f'{folder}_unrelaxed_rank_001.pdb'
                name_new_json = f'{folder}_scores_rank_001.json'

                # copy pdb file of best-ranking model
                path_new_pdb = f'{specification_output}/{self.prediction_folder_name}/{name_new_pdb}'
                shutil.copy(best_pdb, path_new_pdb)

                # produce new score file (in json format) 
                path_new_json = f'{specification_output}/{self.prediction_folder_name}/{name_new_json}'
                format_scores(best_pkl, path_new_json)