import os
import glob
import cPickle
from rdkit import Chem
from rdkit.Chem import AllChem


os.chdir("200_hldiff_1")
gen_list = glob.glob("*_gen.p")
results_list = glob.glob("*_results.p")

uniqs = []
for g in gen_list:
    mol_list = cPickle.load( open(g))
    for mol in mol_list:
        if mol.GetProp("_Energy") == "None":
             uniqs.append(mol)
for r in results_list:
    mol_list = cPickle.load( open(r))
    for mol in mol_list:
        if float(mol.GetProp("_Energy")) > -0.06690:
            print r

print len(uniqs) 
