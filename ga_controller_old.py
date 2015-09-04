import cPickle
import argparse
import ga_crossover
import ga_setup_selection
import glob
from rdkit import Chem
from rdkit.Chem import AllChem

parser = argparse.ArgumentParser(description="-"*78 + "Genetic algorithm controller".ljust(78) + "-"*78)
parser.add_argument("-pop", "--population_size", type=int, help="Size of initial population", default=100)
parser.add_argument("-n_number", "--nitrogen_number", type=int, help="Number of intial nitrogens placed into molecule", default=5)
parser.add_argument("-x_type", "--crossover_type", type=str, help="Type of crossover used, either fixed N or varying", default="varying")
parser.add_argument("-x_rate", "--crossover_rate", type=int, help="amount of each generation produced by crossover (%%)", default=90)
parser.add_argument("-e_rate", "--elitism_rate", type=int, help="amount of each generation produced by elitism (%%)", default=10)
#parser.add_argument("-setup", "--setup_ga_run", type=str, help="setup new run or continue current one", default="Y")
parser.add_argument("-n_gen", "--number_of_generations", type=int, help="how many generations to run for", default=5)

args = parser.parse_args()



gen_counter = 1
gen_name = str(args.population_size)+"_"+str(gen_counter)
smiles_string = 'c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1'
#temp_list = list(smiles_string)
#mutant = ga_setup_selection.mutator(temp_list, mutation_factor = 10, N_number=5)
#mutant_string = ''.join(mutant)
#m3 = Chem.MolFromSmiles(mutant_string)
#print mutant_string
#if args.setup_ga_run == "Y":
    pop = ga_setup_selection.generate_initial_population(smiles_string,args.population_size,args.nitrogen_number)
    print len(pop)
    mol_lst = ga_setup_selection.mol_name_maker(pop)
    all_mol_coords = ga_setup_selection.coords_generator(mol_lst)
    print len(all_mol_coords)
    ga_setup_selection.com_generator_neutral_grndst(all_mol_coords)
    #ga_setup_selection.com_generator_neg_mol_neg_geom(all_mol_coords)
    dump_list = []
    for m in mol_lst:
        name = m.GetProp("_Name")
        pm = AllChem.PropertyMol(m)
        pm.SetProp("_Name", name)
        pm.SetProp("_Energy", None)
        dump_list.append(pm)
    cPickle.dump(dump_list, open(gen_name)+"_gen.p", "w+")

#if args.setup_ga_run == "N":
   #getting right pickle, should always be the latest one
   gen_picks = glob.glob("*_gen.p")
   sorted_gen = sorted(gen_picks, key=lambda gen: int(gen.split("_")[1]))
   latest_gen = sorted_gen[-1]
   pop_size = latest_gen.split("_")[0]
   gen_counter = latest_gen.split("_")[1]
   #opening up latest pickle and results file
   mol_name_list = []
   name_list = cPickle.load( open(latest_gen, "rb" ) )
   for mol in name_list:
       temp_m = mol
       #print temp_m.GetProp("_Name")
       #print temp_m.GetProp("_Energy")
       mol_name_list.append(temp_m)
   f = open("reorg_results_"+gen_counter+"_"+pop_size+".txt", "r")
   res = []
   for line in f:
       lin = line.strip().split()
       namo= lin[0].split("_")[:2]
       name  = "_".join(namo)
       en = lin[1]
       res.append((name,en)) 
   print len(res)
   #matching up none energy mols with results
   mol_en = []
   i = 0
   for mol in mol_name_list:
       if mol.GetProp("_Energy") == "None":
           #print mol.GetProp("_Name"), i
           try:
               mol.SetProp("_Energy",res[i][1])
               i += 1
           except IndexError:
               print "whoops"
       mol_en.append((mol,mol.GetProp("_Energy")))
   print len(mol_en)
   ml_sorted_by_en = sorted(mol_en, key=lambda tup: float(tup[1]))
   #writing results pickle
   dump_list = []
   for mol,en in ml_sorted_by_en:
       name = mol.GetProp("_Name")
       en = mol.GetProp("_Energy")
       pm = AllChem.PropertyMol(mol)
       pm.SetProp("_Name", name)
       pm.SetProp("_Energy", en)
       dump_list.append(pm)
   cPickle.dump(dump_list, open(pop_size+"_"+gen_counter+"_results.p", "w+"))
   #setting crossover/elitism sizes
   crossover_size = int(pop_size) * args.crossover_rate/100
   elite_size = int(pop_size) * args.elitism_rate/100
   print "Using crossover size of "+str(crossover_size)+" and elitism size of "+str(elite_size)
   #generating the right lists, mols are picked for crossover like first with last second with second last etc
   crossover_mols,elite_mols = ga_setup_selection.tournament_selection(mol_list_en = ml_sorted_by_en, elitism_rate=elite_size,new_pop_size=crossover_size)
   pairs_set = []
   pairs = zip(crossover_mols[::1], crossover_mols[-1::-1])
   for p in pairs:
       if (p[0],p[1]) in pairs_set:
           pass
       if (p[1],p[0]) in pairs_set:
           pass
       else:
           pairs_set.append((p[0],p[1]))
   print len(pairs_set)
   children = []
   for p in pairs_set:
       kids = ga_crossover.crossover_pair_varying(p[0],p[1])
       children+=kids
   print len(children)
   #combining elite and crossed over mols
   mol_list_2 = []
   new_gen = elite_mols + children
   for i,ng in enumerate(new_gen):
       ng.SetProp("_Name", str(i)+"_test")
       mol_list_2.append(ng)
   print len(mol_list_2)
   #removing duplicates and writing new generation
   dump_list = ga.setup_selection.check_previous_gens(mol_list_2)
   for m in new_matches_sorted:
    name = m.GetProp("_Name")
    try:
        energy = m.GetProp("_Energy")
    except KeyError:
        en = "None"
        energy = m.SetProp("_Energy", en)
    pm = AllChem.PropertyMol(m)
    pm.SetProp("_Name", name)
    pm.SetProp("_Energy", energy)
    dump_list.append(pm)
    cPickle.dump(dump_list, open("2_9_gen_50.p", "w+"))


