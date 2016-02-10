import cPickle
import argparse
import os
import glob
from rdkit import Chem
from rdkit.Chem import AllChem
import celery
from celery import group
from tasks import run_gauss, run_gauss2
from celery.task.control import broadcast
from celery.task.control import discard_all
import ga_crossover
import ga_setup_selection
import ga_ring_crossover
import ga_gauss

parser = argparse.ArgumentParser(description="-"*78 + "Genetic algorithm controller".ljust(78) + "-"*78)
parser.add_argument("-pop", "--population_size", type=int, help="Size of initial population", default=100)
parser.add_argument("-n_number", "--nitrogen_number", type=int, help="Number of intial nitrogens placed into molecule", default=5)
parser.add_argument("-x_type", "--crossover_type", type=str, help="Type of crossover used, either fixed N/varying or ring crossover", default="ring")
parser.add_argument("-x_rate", "--crossover_rate", type=int, help="amount of each generation produced by crossover (%%)", default=90)
parser.add_argument("-e_rate", "--elitism_rate", type=int, help="amount of each generation produced by elitism (%%)", default=10)
parser.add_argument("-n_gen", "--number_of_gens", type=int, help="how many generations to run for", default=50)
parser.add_argument("-restart", "--restart_ga", type=int, help="restart ga from this generation",default=None)
parser.add_argument("-gauss_calc", "--gauss_calc",type=str,help="which molecular property you want to optimise, reorg_en, hl_diff, max_dipole,lumo",default="reorg_en")
parser.add_argument("-work_dir", "--work_dir",type=str,help="where to run the GA",default=".")


args = parser.parse_args()
number_of_gens = range(1,args.number_of_gens+1)
if args.restart_ga:
   number_of_gens = range(args.restart_ga, args.number_of_gens+1)
print "Beginning ga run with "+str(args.population_size)+" pop members for "+str(args.number_of_gens)+" generations"
try:
    os.chdir(args.work_dir)
except OSError:
    os.mkdir(args.work_dir)
    os.chdir(args.work_dir)
for ng in number_of_gens:
    gen_name = str(args.population_size)+"_"+str(ng)
    print "Running "+gen_name
    p_name = gen_name+"_gen.p"
    if ng == 1:
        #initial setup, pop_size, mol_coords and com_files
        smiles_string = 'c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1'
        pop = ga_setup_selection.generate_initial_population(smiles_string,args.population_size,args.nitrogen_number)
        #print len(pop)
        mol_lst = ga_setup_selection.mol_name_maker(pop)
        #writing generation .p
        ga_setup_selection.init_gen_dumper(p_name,mol_lst)
        #collecting and running coms
        all_mol_coords = ga_setup_selection.coords_generator(mol_lst)
        if args.gauss_calc == "reorg_en":
            gen_result = ga_gauss.reorg_en_calc(all_mol_coords)
        if args.gauss_calc == "hl_diff":
            gen_result = ga_gauss.hl_diff_calc(all_mol_coords)
        if args.gauss_calc == "max_dipole":
            gen_result = ga_gauss.max_dipole_calc(all_mol_coords)
        if args.gauss_calc == "lumo":
            gen_result = ga_gauss.lumo_calc(all_mol_coords)
        print gen_result
        #matching up mol objects and energies
        mol_name_list = ga_setup_selection.gen_loader(p_name)
        if args.gauss_calc == "lumo" or "max_dipole" or "hl_diff":
            ml_sorted_by_en = ga_setup_selection.result_matcher(gen_result,mol_name_list,maximise="True")
        if args.gauss_calc == "reorg_en":
            ml_sorted_by_en = ga_setup_selection.result_matcher(gen_result,mol_name_list,maximise="False") 
        print ml_sorted_by_en
        #writing results pickle
        ga_setup_selection.results_dumper(p_name,ml_sorted_by_en)
        #setting crossover/elitism sizes
        crossover_size = int(args.population_size) * args.crossover_rate/100
        elite_size = int(args.population_size) * args.elitism_rate/100
        print "Using crossover size of "+str(crossover_size)+" and elitism size of "+str(elite_size)
        #generating the right lists, mols are picked for crossover like first with last second with second last etc
        crossover_mols,elite_mols = ga_setup_selection.tournament_selection(mol_list_en = ml_sorted_by_en, elitism_rate=elite_size,new_pop_size=crossover_size)
        pairs_set = ga_setup_selection.make_pairs(crossover_mols)
        children = []
        for p in pairs_set:
            while True:
                try:
                    kids = ga_ring_crossover.ring_crossover(p[0],p[1])
                    children+=kids
                    break
                except (IndexError,RuntimeError,ValueError):
                    continue  
        #combining elite and crossed over mols
        mol_list_2 = []
        new_gen = elite_mols + children
        for i,nm in enumerate(new_gen,start=1):
            nm.SetProp("_Name", str(i)+"_test")
            mol_list_2.append(nm)
        #removing duplicates and writing new generation
        new_mols = ga_setup_selection.check_previous_gens(mol_list_2)
        ga_setup_selection.new_gen_dumper(args.population_size,ng,new_mols)
    else:
        mol_list = cPickle.load(open(p_name, "rb" ))
        print len(mol_list)
        run_mols = [mol for mol in mol_list if mol.GetProp("_Energy") == "None"]
        print len(run_mols)
        if len(run_mols) <= 1:
            print "GA has hit generation of size 1, will quit now"
            print "Shutting down workers"
            broadcast("shutdown")
            exit()
        all_mol_coords = ga_setup_selection.coords_generator(run_mols)
        if args.gauss_calc == "reorg_en":
            gen_result = ga_gauss.reorg_en_calc(all_mol_coords)
        if args.gauss_calc == "hl_diff":
            gen_result = ga_gauss.hl_diff_calc(all_mol_coords)
        if args.gauss_calc == "max_dipole":
            gen_result = ga_gauss.max_dipole_calc(all_mol_coords)
        if args.gauss_calc == "lumo":
            gen_result = ga_gauss.lumo_calc(all_mol_coords)
        print gen_result
        if args.gauss_calc == "lumo" or "max_dipole" or "hl_diff":
            ml_sorted_by_en = ga_setup_selection.result_matcher(gen_result,mol_list,maximise="True")
        if args.gauss_calc == "reorg_en":
            ml_sorted_by_en = ga_setup_selection.result_matcher(gen_result,mol_list,maximise="False")
        print ml_sorted_by_en
        #writing results pickle
        ga_setup_selection.results_dumper(p_name,ml_sorted_by_en)
        new_min = float(ml_sorted_by_en[0][1])
        if new_min < 0.1333:
            print "reached min, stopping"
            broadcast("shutdown")
            exit()
        #setting crossover/elitism sizes
        crossover_size = int(args.population_size) * args.crossover_rate/100
        elite_size = int(args.population_size) * args.elitism_rate/100
        print "Using crossover size of "+str(crossover_size)+" and elitism size of "+str(elite_size)
        #generating the right lists, mols are picked for crossover like first with last second with second last etc
        crossover_mols,elite_mols = ga_setup_selection.tournament_selection(mol_list_en = ml_sorted_by_en, elitism_rate=elite_size,new_pop_size=crossover_size)
        pairs_set = ga_setup_selection.make_pairs(crossover_mols)
        children = []
        for p in pairs_set:
            while True:
                try:
                    kids = ga_ring_crossover.ring_crossover(p[0],p[1])
                    children+=kids
                    break
                except (IndexError,RuntimeError,ValueError):
                    continue
        print len(children)
        #combining elite and crossed over mols
        mol_list_2 = []
        new_gen = elite_mols + children
        for i,nm in enumerate(new_gen,start=1):
            nm.SetProp("_Name", str(i)+"_test")
            mol_list_2.append(nm)
        #removing duplicates and writing new generation
        new_mols = ga_setup_selection.check_previous_gens(mol_list_2)
        ga_setup_selection.new_gen_dumper(args.population_size,ng,new_mols)

