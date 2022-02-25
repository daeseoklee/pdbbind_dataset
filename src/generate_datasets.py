from pathlib import Path
from this import d 
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from tqdm import tqdm
import json, pickle
import shutil
import numpy as np


"""
This module is supposed to be run after running
    -initial parts of pdbbind_data.ipynb 
    -split_data.py 
    -clean_pdbs.ipynb
"""

class SDFParsingException(Exception):
    def __init__(self, sdf_file):
        super().__init__()
        self.sdf_file = sdf_file
    def __str__(self):
        return f'Error while parsing {self.sdf_file}'

class SMILESGenerationError(Exception):
    def __init__(self, sdf_file):
        super().__init__()
        self.sdf_file = sdf_file
    def __str__(self):
        return f'Error while obtaining rdkit-canonical SMILES molecule in {self.sdf_file}'

def get_canonical_permutation(sdf_file):
    count = 0
    for mol in Chem.SDMolSupplier(sdf_file, sanitize=False):
        count += 1
    assert count == 1
    if mol is None:
        raise SDFParsingException(sdf_file)
    mol.UpdatePropertyCache(strict=False)
    mol = Chem.rdmolops.RemoveAllHs(mol, sanitize=False) 
    try:
        target_smiles = Chem.MolToSmiles(mol)
    except:
        raise SMILESGenerationError(sdf_file) 
    target_mol = Chem.MolFromSmiles(target_smiles, sanitize=False) #will reorder as this 
    #try:
    #    Chem.rdmolops.SanitizeMol(target_mol)
    #except:
    #    return target_smiles, None
    if target_mol is None:
        raise Exception()
    target_mol.UpdatePropertyCache(strict=False)
    
    orig_ranks = list(Chem.CanonicalRankAtoms(mol))
    orig_ranks_reversed = tuple(zip(*sorted([(j, i) for i, j in enumerate(orig_ranks)])))[1]

    target_ranks = Chem.CanonicalRankAtoms(target_mol)
    if not len(orig_ranks_reversed) == len(target_ranks):
        raise Exception()
    perm = [orig_ranks_reversed[rank] for rank in target_ranks]
    mol_renum = Chem.RenumberAtoms(mol, perm)
    for atom1, atom2 in zip(mol_renum.GetAtoms(), target_mol.GetAtoms()):
        if atom1.GetSymbol() != atom2.GetSymbol():
            raise Exception()
    return target_smiles, perm

def get_atom_coord_list(sdf_file, perm):
    orig_atom_coord_list = []
    def get_begin_idx(lines):
        return 4
    def get_end_idx(lines):
        for i, line in enumerate(lines[4:], start=4):
            if not '.' in line:
                return i - 1
    with open(sdf_file, 'r') as f:
        lines = f.readlines()
        begin_idx = get_begin_idx(lines)
        end_idx = get_end_idx(lines)
        for line in lines[begin_idx:end_idx+1]:
            strs = line.split()
            x, y, z, symbol, *_ = strs 
            x, y, z = float(x), float(y), float(z)
            if symbol == 'H':
                continue
            orig_atom_coord_list.append((symbol, (x, y, z)))
    if len(perm) != len(orig_atom_coord_list):
        msg = '--------------------\n'
        msg += f'len(perm): {len(perm)}\n'
        msg += f'len(orig_atom_coord_list): {len(orig_atom_coord_list)}\n'
        msg += f'See: {sdf_file}\n'
        msg += f'begin_idx: {begin_idx}\n'
        msg += f'end_idx: {end_idx}\n'
        for atomname, coord in orig_atom_coord_list:
            msg += f'{atomname} - {coord}\n'
        raise Exception(msg)
    atom_coord_list = [None for _ in range(len(perm))]
    for i in range(len(perm)):
        atom_coord_list[i] = orig_atom_coord_list[perm[i]]
    return atom_coord_list

def check_aligned(atom_coord_list, target_smiles):
    target_mol = Chem.MolFromSmiles(target_smiles, sanitize=False)
    for (atomname, coord), atom in zip(atom_coord_list, target_mol.GetAtoms()):
        assert atomname == atom.GetSymbol()

def get_fasta_list(data_dir_list):
    print('Obtaining fasta..')

    with open(cwd / 'cc_fasta_dict.json', 'r') as f:
        cc_fasta_dict = json.load(f)
        
    fasta_list = []
    for data_dir in tqdm(data_dir_list):
        pdb_code = data_dir.name
        cc_fasta = cc_fasta_dict[pdb_code]
        fasta_list.append(cc_fasta)
    return fasta_list

class SanitizationError(Exception):
    def __init__(self):
        super().__init__()

def get_smiles_and_atom_coord_list(data_dir):
    sdf_file = str(data_dir / f'{data_dir.name}_ligand.sdf') 
    target_smiles, perm = get_canonical_permutation(sdf_file)
    
    if perm is None:
        raise SanitizationError()
    
    atom_coord_list = get_atom_coord_list(sdf_file, perm)
    check_aligned(atom_coord_list, target_smiles)
    
    return target_smiles, atom_coord_list
        
def get_smiles_list(data_dir_list):
    print('Obtaining smiles..')
    smiles_list = []
    sdf_parsing_exceptions = 0
    smiles_exceptions = 0
    sanitization_errors  = 0
    for data_dir in tqdm(data_dir_list):
        try:
            target_smiles, _ = get_smiles_and_atom_coord_list(data_dir)
        except SDFParsingException:
            sdf_parsing_exceptions += 1
            smiles_list.append(None)
            #raise Exception(sdf_file)
            raise Exception(data_dir)
            continue
        except SMILESGenerationError:
            smiles_exceptions += 1
            smiles_list.append(None)
            continue
        except SanitizationError:
            sanitization_errors += 1
            smiles_list.append(None)
            continue
        except:
            raise Exception()
        smiles_list.append(target_smiles)
    print(f'sdf parsing exceptions: {sdf_parsing_exceptions}')
    print(f'SMILES generation errors: {smiles_exceptions}')
    print(f'sanitization errors: {sanitization_errors}')
    
    return smiles_list
    
def _gather_dataset_info(split):
    def filter_lists(fasta_list, smiles_list, affinity_list, data_dir_list):
        """
        Smiles == None is excluded 
        """
        new_fasta_list = []
        new_smiles_list = []
        new_affinity_list = []
        new_data_dir_list = []
        for (fasta, smiles, affinity, data_dir) in zip(fasta_list, smiles_list, affinity_list, data_dir_list):
            if smiles is None:
                continue 
            new_fasta_list.append(fasta)
            new_smiles_list.append(smiles)
            new_affinity_list.append(affinity)
            new_data_dir_list.append(data_dir)
        return new_fasta_list, new_smiles_list, new_affinity_list, new_data_dir_list
            
    with open(dataset_dir / f'{split}_data.tsv', 'r') as f:
        affinity_list = [] 
        data_dir_list = []
        for line in tqdm(f.readlines()):
            code, affinity, pdbbind_set = line.split()
            affinity_list.append(affinity)
            data_dir_list.append(pdbbind_set_to_dir[pdbbind_set] / code)
    fasta_list = get_fasta_list(data_dir_list)
    #fasta_list = None
    smiles_list = get_smiles_list(data_dir_list)
    
    fasta_list, smiles_list, affinity_list, data_dir_list = filter_lists(fasta_list, smiles_list, affinity_list, data_dir_list)
    assert len(fasta_list) == len(smiles_list) == len(affinity_list) == len(data_dir_list)
    
    fasta_set = set(fasta_list)
    smiles_set = set(smiles_list)
    triple_list = list(zip(fasta_list, smiles_list, affinity_list))
            
    return fasta_set, smiles_set, triple_list, data_dir_list


def parse_clean_pdb(clean_pdb, return_coords=True):

    coords = []
    ca_exist = [] 
    with open(clean_pdb, 'r') as reader:
        
        chainID, resSeq, iCode = None, None, None
        CA_encountered = True
        for line in reader:
            assert line.startswith('ATOM') or line.startswith('HETATM')
            prev_chainID = chainID 
            prev_resSeq = resSeq 
            prev_iCode = iCode
            
            chainID = line[21]
            resSeq = int(line[22:26])
            iCode = line[26]
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            atomname = line[12:16].strip()
                    
            
            at_new_residue = (chainID, resSeq, iCode) != (prev_chainID, prev_resSeq, prev_iCode)
            
            if at_new_residue:
                if coords != []:
                    ca_exist.append(CA_encountered)
                CA_encountered = False
                coords.append([[x, y, z]])
            else:
                coords[-1].append([x, y, z]) 
            
            if atomname == 'CA':
                CA_encountered = True
            
        ca_exist.append(CA_encountered)
    
    ca_exist = np.array(ca_exist)
    
    if not return_coords:
        return ca_exist
    
    max_len = max(len(res_coords) for res_coords in coords)
    coords_mask = [[True] * len(res_coords) + [False] * (max_len - len(res_coords)) for res_coords in coords]
    coords = [res_coords + [[0.,0.,0.] for _ in range(max_len - len(res_coords))] for res_coords in coords]
    coords_mask = np.array(coords_mask)
    coords = np.array(coords)
    
    return ca_exist, coords, coords_mask
    

def generate_pdbbind_ccsan_dti_datasets():
    split_to_fasta_set = {}
    split_to_smiles_set = {}
    split_to_triple_list = {}
    split_to_data_dir_list = {}
    for split in ['train', 'dev', 'test']: 
        print(f'<Preprocessing {split}>')
        split_to_fasta_set[split], split_to_smiles_set[split], split_to_triple_list[split], split_to_data_dir_list[split] = _gather_dataset_info(split)

    fasta_list = list(split_to_fasta_set['train'] | split_to_fasta_set['dev'] | split_to_fasta_set['test'])
    smiles_list = list(split_to_smiles_set['train'] | split_to_smiles_set['dev'] | split_to_smiles_set['test'])
    
    with open(pdbbind_ccsan_derived_dir / 'ligand_index.json', 'w') as f:
        ligand_index_dict = {
            str(i): smiles for i, smiles in enumerate(smiles_list)
        }
        json.dump(ligand_index_dict, f)
    
    with open(pdbbind_ccsan_derived_dir / 'protein_index.json', 'w') as f:
        protein_index_dict = {
            str(i): fasta for i, fasta in enumerate(fasta_list)
        }
        json.dump(protein_index_dict, f)
    
    for split in ['train', 'dev', 'test']: 
        print(f'<processing {split}>')
        _triple_list = split_to_triple_list[split]
        triple_list = [(str(fasta_list.index(fasta)), str(smiles_list.index(smiles)), affinity) for fasta, smiles, affinity in _triple_list]
        with open(pdbbind_ccsan_derived_dir / f'pdbbind_ccsan_pafnucy_{split}.tsv', 'w') as f:
            for p_idx, m_idx, affinity in triple_list:
                f.write(f'{p_idx}\t{m_idx}\t{affinity}\n')
        
        data_dir_list = split_to_data_dir_list[split]
        
        with open(pdbbind_ccsan_derived_dir / f'pafnucy_{split}_pdbcode_list.pkl', 'wb') as f:
            code_list = [data_dir.name for data_dir in data_dir_list]
            pickle.dump(code_list, f)
        
        clean_pdb_target_dir = pdbbind_ccsan_derived_dir / 'clean_pdbs' 
        if not clean_pdb_target_dir.exists():
            clean_pdb_target_dir.mkdir()
        pdbcode_to_resmask_acc = {}
        resmask_acc_list = []
        resmask_list = []
        compact_dist_matrix_list = [] 
        dist_matrix_list = [] 
        print('Copying and precomputing for clean pdbs...')
        for data_dir in tqdm(data_dir_list):
            clean_pdb_file = data_dir / f'{data_dir.name}_clean.pdb'
            clean_pdb_target_file = clean_pdb_target_dir / f'{data_dir.name}.pdb'
            if not clean_pdb_target_file.exists():
                shutil.copyfile(str(clean_pdb_file), str(clean_pdb_target_file))
            
            resmask = parse_clean_pdb(clean_pdb_target_file, return_coords=False)
            resmask_acc = (np.tril(np.ones((len(resmask), len(resmask)))) @ resmask).astype(np.int) 
            
            pdbcode_to_resmask_acc[data_dir.name] = resmask_acc 
            resmask_acc_list.append(resmask_acc)
            resmask_list.append(resmask)
        
            #for compact_dist_matrix_list

            pdbcode = data_dir.name
            clean_pdb_file = data_dir / f'{pdbcode}_clean.pdb'
            
            try:
                _, atom_coord_list = get_smiles_and_atom_coord_list(data_dir)
            except:
                raise Exception('could not parse ligand coords:', data_dir)

            ca_exist, coords, coords_mask = parse_clean_pdb(clean_pdb_file, return_coords=True)
            
            ligand_coords = np.array([[x, y, z] for _, (x, y, z) in atom_coord_list])
            
            diff = coords[:, :, None, :] - ligand_coords[None, None, :, :]
            dist = np.sqrt(np.sum(diff ** 2, axis=-1))
            dist_matrix = np.min(np.where(coords_mask[:, :, None], dist, np.inf), axis=1) 
            #residue-atom distances, (if no residue coordinate is provided, infinity)
            
            tentimes_dist_matrix = 10 * dist_matrix
            dist_matrix = np.uint8(np.where(tentimes_dist_matrix<255, tentimes_dist_matrix, 255))
            
            compact_dist_matrix = dist_matrix[ca_exist, :]
                
            compact_dist_matrix_list.append(compact_dist_matrix)
            dist_matrix_list.append(dist_matrix_list)
            
        with open(pdbbind_ccsan_derived_dir / f'pafnucy_{split}_resmask_acc_list.pkl', 'wb') as f:
            pickle.dump(resmask_acc_list, f)

        with open(pdbbind_ccsan_derived_dir / f'pafnucy_{split}_resmask_list.pkl', 'wb') as f:
            pickle.dump(resmask_list, f)
        
        with open(pdbbind_ccsan_derived_dir / f'pafnucy_{split}_compact_dist_matrix_list.pkl', 'wb') as f:
            pickle.dump(compact_dist_matrix_list, f)
            
        with open(pdbbind_ccsan_derived_dir / f'pafnucy_{split}_dist_matrix_list.pkl', 'wb') as f:
            pickle.dump(dist_matrix_list, f)
            

def _generate_pbsdb_int_dataset(split, collect_clean_pdbs=True):
    """
    1. Collect (data_dir) candidates from dataset_dir/{split}.tsv
    2. Try obtaining canonical SMILES as well as coordinates using "get_smiles_and_atom_coord_list(data_dir)" 
    3. If fail, goto step2 for the next (code, set)
    4. copy the clean pdb to pbsdb/clean_pdbs/ 
    5. Add (fasta, smiles) |-> (intdist_dict, None, chainstr, None, res_mask) to a dict, where 
        -intdist_dict: 
        -chainstr: the pdb code
        -res_mask: mask of positions where CA is provided 
        -None: placeholder to match the corresponding information in kiba
        
        When (fasta, smiles) already exists, skip the instance and count such event
    6. After finishing the loop, pickle the dict to pbsdb/{split}_data_dict.pkl 
    """
    print(f'<pbsdb {split}>')
    def loop_over_data_dir():
        with open(dataset_dir / f'{split}_data.tsv', 'r') as f:
            for line in f:
                code, _, pdbbind_set = line.split()
                data_dir = pdbbind_set_to_dir[pdbbind_set] / code
                yield data_dir

    clean_pdb_target_dir = pbsdb_dir / 'clean_pdbs'
    if not clean_pdb_target_dir.exists():
        clean_pdb_target_dir.mkdir()

    with open(cwd / 'cc_fasta_dict.json', 'r') as f:
        cc_fasta_dict = json.load(f)

    sdf_parsing_exceptions = 0
    smiles_errors = 0
    sanitization_errors  = 0
    pair_duplicates = 0 
    data_dict = {}
    
    for data_dir in tqdm(list(loop_over_data_dir())):
        pdb_code = data_dir.name
        fasta = cc_fasta_dict[pdb_code]
        try:
            smiles, atom_coord_list = get_smiles_and_atom_coord_list(data_dir)
        except SDFParsingException:
            sdf_parsing_exceptions += 1
            continue 
        except SMILESGenerationError:
            smiles_errors += 1
        except SanitizationError:
            sanitization_errors += 1
            continue
        except:
            raise Exception()
        if (fasta, smiles) in data_dict:
            pair_duplicates += 1
            continue 
        
        clean_pdb_file = data_dir / f'{pdb_code}_clean.pdb'
        clean_pdb_file_target =  clean_pdb_target_dir / f'{pdb_code}.pdb'
        shutil.copyfile(str(clean_pdb_file), str(clean_pdb_file_target))
        
        
        ca_exist, coords, coords_mask = parse_clean_pdb(clean_pdb_file, return_coords=True)
        assert len(fasta) == len(coords) == len(ca_exist) 
        
        ligand_coords = np.array([[x, y, z] for _, (x, y, z) in atom_coord_list])
        
        diff = coords[:, :, None, :] - ligand_coords[None, None, :, :]
        dist = np.sqrt(np.sum(diff ** 2, axis=-1))
        dist_matrix = np.min(np.where(coords_mask[:, :, None], dist, np.inf), axis=1) 
        #residue-atom distances, (if no residue coordinate is provided, infinity)
        intdist_dict = {}
        for i, j in zip(*np.where(dist_matrix < 7.0)):
            intdist_dict[(i,j)] = dist_matrix[i, j] 
            
        resmask = ca_exist
        
        data_dict[(fasta, smiles)] = (intdist_dict, None, pdb_code, None, resmask)
    
    with open(pbsdb_dir / f'{split}_data_dict.pkl', 'wb') as f:
        pickle.dump(data_dict, f)
    
    print('sdf parsing exceptions:', sdf_parsing_exceptions)
    print('SMILES generation errors:', smiles_errors)
    print('sanitization errors:', sanitization_errors)
    print('pair duplicates', pair_duplicates)
        
        

def generate_pbsdb_int_datasets(collect_clean_pdbs=True):
    for split in ['train', 'dev', 'test']:
        _generate_pbsdb_int_dataset(split, collect_clean_pdbs=collect_clean_pdbs)
    


if __name__ == '__main__':
    cwd = Path(__file__).parent
    dataset_dir = cwd.parent / 'dataset'
    
    pdbbind_base_dir = cwd.parent / 'pdbbind' / 'v2016'
    pdbbind_set_to_dir = {
        'general': pdbbind_base_dir / 'general-set-except-refined',
        'refined': pdbbind_base_dir / 'refined-set',
        'core': pdbbind_base_dir / 'refined-set'
    }
    
    
    
    #dti dataset--------------------------------------------------------------------------
    pdbbind_ccsan_derived_dir = cwd.parent / 'pdbbind_ccsan_derived'
    if not pdbbind_ccsan_derived_dir.exists():
        pdbbind_ccsan_derived_dir.mkdir()
    
    generate_pdbbind_ccsan_dti_datasets()
    
    
    #int dataset --------------------------------------------------------------------------
    pbsdb_dir = cwd.parent / 'pbsdb'
    if not pbsdb_dir.exists():
        pbsdb_dir.mkdir()
    
    generate_pbsdb_int_datasets()    
    
    
    
