import pandas as pd
import csv
from sklearn.utils import shuffle
import numpy as np

def remove_core2013(d):
    code_set_list = list(zip(d['pdbid'], d['set']))
    general = set([code for code, s in code_set_list if s == 'general' or s == 'refined'])
    core = set([code for code, s in code_set_list if s == 'core'])
    
    with open('dataset/core_pdbbind2013.ids', 'r') as f:
        core2013 = set([line.strip() for line in f])
    d['include'] = True
    d.loc[np.in1d(d['pdbid'], list(core2013 & (general - core))), 'include'] = False
    return d[d.include][['pdbid', '-logKd/Ki', 'set']]

def assign_splits_(d):
    d['split'] = 'train' 
    d.loc[d['set'] == 'core', 'split'] = 'test'
    
    in_refined = np.array(d.set == 'refined')
    indices = np.where(in_refined)[0]
    shuffled_indices = shuffle(list(indices), random_state=123)
    sampled_indices = np.array(sorted(shuffled_indices)[:1000])
    in_dev = np.zeros(len(in_refined), dtype=bool) 
    in_dev[sampled_indices] = True
    d.loc[in_dev, 'split'] = 'dev'

def print_stats(d):
    print(d.groupby('split').apply(len).loc[['train', 'dev', 'test']])
    print(d.head())

if __name__ == '__main__':
    d = pd.read_csv('dataset/affinity_data_cleaned.csv')
    d = remove_core2013(d)
    assign_splits_(d)
    print_stats(d)
    
    for split in ['train', 'test', 'dev']:
        d.loc[d.split == split][['pdbid', '-logKd/Ki', 'set']].to_csv(f'dataset/{split}_data.tsv', header=False, sep='\t', index=False)
    
    
    
    


    