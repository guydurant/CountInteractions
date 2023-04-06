from pymol import cmd
import os
from tqdm import tqdm
from plip import structure
from plip.structure import preparation, detection
from plip.structure.preparation import PDBComplex
import pandas as pd

def load_csv(csv_file, data_dir):
    df = pd.read_csv(csv_file)
    protein_files = [os.path.join(data_dir, protein_file) for protein_file in df['protein']]
    ligand_files = [os.path.join(data_dir, ligand_file) for ligand_file in df['ligand']]
    keys = df['key']
    pks = df['pk']
    return protein_files, ligand_files, keys, pks

def combine_protein_ligand_files(protein_file, ligand_file, pdb):
    cmd.load(protein_file, 'protein')
    cmd.load(ligand_file, 'ligand')
    cmd.select('complex', 'all')
    cmd.create('complex_obj', 'complex')
    cmd.save(f'temp_files/{pdb}_complex.pdb', 'complex_obj')
    cmd.delete('all')
    return None

def count_interactions(protein_file, ligand_file, pdb, pk):
    combine_protein_ligand_files(protein_file, ligand_file, pdb)
    my_mol = PDBComplex()
    my_mol.load_pdb(f'temp_files/{pdb}_complex.pdb', as_string=False)
    my_mol.analyze()
    my_interactions = my_mol.interaction_sets
    key = list(my_interactions.keys())[0]
    return {
        'total': len(my_interactions[key].all_itypes),
        'hbonds': len(my_interactions[key].hbonds_ldon)+len(my_interactions[key].hbonds_pdon),
        'hydrophobic': len(my_interactions[key].hydrophobic_contacts),
        'pk': pk,
    }

def analyse_all_interactions(csv_file, data_dir):
    protein_files, ligand_files, keys, pks = load_csv(csv_file, data_dir)
    results = {}
    for i in tqdm(range(len(protein_files))):
        results[keys[i]] = count_interactions(protein_files[i], ligand_files[i], keys[i], pks[i])
        os.remove(f'temp_files/{keys[i]}_complex.pdb')
    return pd.DataFrame(results)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv_file', type=str, default='train.csv')
    parser.add_argument('--val_csv_file', type=str, default='val.csv')
    parser.add_argument('--data_dir', type=str, default='data')
    parser.add_argument('--val_data_dir', type=str, default='data')
    parser.add_argument('--model_name', type=str, default='test')
    parser.add_argument('--train', action='store_true')
    parser.add_argument('--predict', action='store_true')
    args = parser.parse_args()
    if args.train:
        raise NotImplementedError('No need to train a model for this task')
    elif args.predict:
        df = analyse_all_interactions(args.val_csv_file, args.val_data_dir)
        df.to_csv(f'results/{args.val_csv_file.split("/")[-1].split(".")[0]}_count_interactions.csv')
    else:
        raise ValueError('Please specify --train or --predict')


