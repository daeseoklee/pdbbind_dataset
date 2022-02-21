import subprocess as sp 
from pathlib import Path
from tqdm import tqdm 

if __name__ == '__main__':
    cwd = Path(__file__).parent
    
    pdbbind_base_dir = cwd.parent / '/home/daelee/pafnucy/pdbbind/v2016'
    pdbbind_dirs = [pdbbind_base_dir / 'general-set-except-refined', pdbbind_base_dir / 'refined-set']
    
    for pdbbind_dir in pdbbind_dirs:
        for data_dir in tqdm(list(pdbbind_dir.glob('*'))):
            data_dir:Path
            if data_dir.name in ['index', 'readme'] or not data_dir.is_dir():
                continue 
            mol2file = str(data_dir / f'{data_dir.name}_ligand.mol2')
            outsmifile = str(data_dir / f'{data_dir.name}_ligand.smi')
            sp.call(['obabel', mol2file, '-O', outsmifile, '-d'])
            #outsdffile = str(data_dir / f'{data_dir.name}_ligand_obabel.sdf')
            #sp.call(['obabel', mol2file, '-O', outsdffile, '-d'])
