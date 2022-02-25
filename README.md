# PDBBind dataset generation 

Files `environment_cpu.yml` and `pdbbind_data.ipynb` are partly from `https://gitlab.com/cheminfIBB/pafnucy`.

## 

## How to generate dataset

1. Prepare PDBBind database by downloading the 2016 version's general set and refined set as `/pdbbind/v2016/general-set-except-refined` and `/pdbbind/v2016/refined-set`.

2. Run
```
python src/split_data.py
```

2. Run the content of `src/clean_pdbs.ipynb`

3. Run 
```
python src/generate_datasets.py 
```

