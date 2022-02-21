# PDBBind dataset generation 

## 

## How to generate dataset
1. Prepare the environment:
```
conda env create -f environment_cpu.yml 
source activate pafnucy_env
```
2. Prepare PDBBind database by downloading the 2016 version's general set and refined set under `/pdbbind/v2016`. (`/pdbbind/v2016/general-set-except-refined` and `/pdbbind/v2016/refined-set`)
3. Run  
```
python src/run_obabel.py #You may have to uncomment some lines 
```
4. Run the content of `src/clean_pdbs.ipynb`

5. Run 
```
python generate_datasets.py #You may have to uncomment some lines
```

