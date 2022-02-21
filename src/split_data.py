import random 
from pathlib import Path

def get_core2013_code_list():
    cwd = Path(__file__).parent
    pdbbind_database_dir = cwd.parent / 'pdbbind' / 'v2016' / 'PDBbind_2016_plain_text_index' / 'index'
    core2013_index_file = pdbbind_database_dir / 'INDEX_core_data.2013'
    l = []
    with open(core2013_index_file, 'r') as f:
        for line in f:
            if line.strip() == '' or line.startswith('#'):
                continue 
            l.append(line.split(' ', 1)[0])
    return l
    

if __name__ == '__main__':
    cwd = Path(__file__).parent
    dataset_dir = cwd.parent / 'dataset'
    data_cleaned_file = dataset_dir / 'affinity_data_cleaned.csv'
    
    all_data = set([])
    with open(data_cleaned_file, 'r') as f:
        f.readline()
        for line in f:
            code, affinity, pdbbind_set = line.strip().split(',')
            all_data.add((code, affinity, pdbbind_set))
    
    orig_d = {'general': set([]), 'refined': set([]), 'core': set([])}
    core2013_code_set = set(get_core2013_code_list())
    num_excluded_general = 0 
    num_excluded_refined = 0
    
    for code, affinity, pdbbind_set in all_data:
        if pdbbind_set in ['general', 'refined']:
            if code in core2013_code_set:
                if pdbbind_set == 'general':
                    num_excluded_general += 1
                if pdbbind_set == 'refined':
                    num_excluded_refined += 1 
                continue
        orig_d[pdbbind_set].add((code, affinity, pdbbind_set)) 
    print('excluded from general\\refined:', num_excluded_general)
    print('excluded from refined:', num_excluded_refined)
    
    general_set, refined_set, core_set = orig_d['general'], orig_d['refined'], orig_d['core']
    dev_set = set(random.sample(list(refined_set), 1000)) 
    train_set = general_set | refined_set.difference(dev_set)
    test_set = core_set 
    split_d = {'train': train_set, 'dev': dev_set, 'test': test_set}
    
    assert num_excluded_general + num_excluded_refined + len(train_set) + len(dev_set) + len(test_set) == len(all_data)
    
    for split, split_set in split_d.items():
        with open(dataset_dir / f'{split}_data.tsv', 'w') as f:
            for code, affinity, pdbbind_set in split_set:
                f.write(f'{code}\t{affinity}\t{pdbbind_set}\n')

    

    