
# MacFrag
MacFrag is an efficient molecule fragmentation method, which is capable of segmenting large-scale molecules in a rapid speed and 
generating diverse fragments that are more compliant with the “Rule of Three”. 

Reference: Yanyan Diao, Feng Hu, Zihao Shen, Honglin Li*. MacFrag: segmenting large-scale molecules to obtain diverse fragments with high qualities. Bioinformatics, 2023. 39(1) : btad012

MacFrag is developed and maintained by Prof. HongLin Li's Group, School of Pharmacy, East China University of Science & Technology, Shanghai 200237, China. 
http://lilab-ecust.cn/

## Modifications

With my current modifications you can also use a mol2 file with a large molecule that you want to fragment. If you have an issue with the code that can possibly be due to rdkit not being able to load your molecule(s), but you should now get a warning if something went wrong. I also added default parameters and made an examples and _data directory so it looks nicer.

## Usage

For help just run 
```py
python macfrag.py --help 
```

## Other files
``` 
1) chembl28_mw500.smi,  Mols collected from ChEMBL database with molecular weight lower than 500
2) chembl28_mw500-1000.smi,  Mols collected from ChEMBL database with molecular weight ranging from 500 to 1000
After merging the two .smi files, 1 921 745 molecules with molecular weights lower than 1000 will be obtained 
that were used to evaluate the qualities of fragments obtained by the three programs. 
3) time_compare.py,  A python script for calculating the run time of different fragmentation programs
``` 
