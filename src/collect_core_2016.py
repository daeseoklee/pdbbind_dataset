from pathlib import Path


def get_core_pdb_codes():
    casf2016_dir = Path('/home/daelee/pdbbind/tmp/CASF-2016')
    core2016_dir = casf2016_dir / 'coreset'
    
    return [path.name for path in core2016_dir.glob('*')]
    

if __name__ == '__main__':
    cwd = Path(__file__).parent
    pdbbind_data_dir = cwd.parent / 'pdbbind' / 'v2016' / 'PDBbind_2016_plain_text_index' / 'index'
    with open(pdbbind_data_dir / 'INDEX_core.2016', 'w') as f:
        for code in get_core_pdb_codes:
            f.write(f'{code}\n')
    
    
    