import cPickle
import argparse
import ga_crossover
import ga_setup_selection
import ga_ring_crossover
import glob
from rdkit import Chem
from rdkit.Chem import AllChem
import celery
from celery import group
from tasks import run_gauss, run_gauss2
from celery.task.control import broadcast
from celery.task.control import discard_all
discard_all()


parser = argparse.ArgumentParser(description="-"*78 + "Genetic algorithm controller".ljust(78) + "-"*78)
parser.add_argument("-pop", "--population_size", type=int, help="Size of initial population", default=100)
parser.add_argument("-n_number", "--nitrogen_number", type=int, help="Number of intial nitrogens placed into molecule", default=5)
parser.add_argument("-x_type", "--crossover_type", type=str, help="Type of crossover used, either fixed N/varying or ring crossover", default="ring")
parser.add_argument("-x_rate", "--crossover_rate", type=int, help="amount of each generation produced by crossover (%%)", default=90)
parser.add_argument("-e_rate", "--elitism_rate", type=int, help="amount of each generation produced by elitism (%%)", default=10)
parser.add_argument("-n_gen", "--number_of_gens", type=int, help="how many generations to run for", default=50)
parser.add_argument("-restart", "--restart_ga", type=int, help="restart ga from this generation",default=None)
parser.add_argument("-gauss_calc", "--gauss_calc",type=str,help="which molecular property you want to optimise, reorg_en, dipole etc",default="reorg_en")



args = parser.parse_args()
number_of_gens = range(1,args.number_of_gens+1)
if args.restart_ga:
   number_of_gens = range(args.restart_ga, args.number_of_gens+1)
print "Beginning ga run with "+str(args.population_size)+" pop members for "+str(args.number_of_gens)+" generations"
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
        dump_list = []
        for m in mol_lst:
            name = m.GetProp("_Name")
            pm = AllChem.PropertyMol(m)
            pm.SetProp("_Name", name)
            pm.SetProp("_Energy", None)
            dump_list.append(pm)
        cPickle.dump(dump_list, open(p_name, "w+"))
        #collecting and running coms
        all_mol_coords = ga_setup_selection.coords_generator(mol_lst)
        com_files_neut = ga_setup_selection.com_generator_neutral_grndst(all_mol_coords)
        com_files_neg = ga_setup_selection.com_generator_neg_mol_neg_geom(all_mol_coords)
        com_files = com_files_neut + com_files_neg
        job = group([run_gauss.s(com) for com in com_files])
        result = job.apply_async()
        result1 = result.get()
        opt_structs = [(name,struct) for name,en,struct in result1 if (name,struct)]
        first_energies = [(name,en) for name,en,struct in result1 if (name,en)]
        #print first_energies
        neg_structs = [(name,struct) for name,struct in opt_structs if "neg" in name]
        neut_structs = [(name,struct) for name,struct in opt_structs if "neg" not in name]
        com_files_neg_2 = ga_setup_selection.com_generator_neg_mol_neut_charge(neg_structs)
        com_files_neut_2 = ga_setup_selection.com_generator_neut_mol_neg_charge(neut_structs) 
        com_files_2 = com_files_neg_2 + com_files_neut_2
        job2 = group([run_gauss2.s(com) for com in com_files_2])
        result2 = job2.apply_async()
        second_energies = result2.get()
        #print first_energies
        #print second_energies
        fe_sorted = sorted(first_energies, key=lambda tup: int(tup[0].split("_")[0]))
        se_sorted = sorted(second_energies, key=lambda tup: int(tup[0].split("_")[0]))
        energies = zip(fe_sorted,se_sorted)
        gen_result = ga_setup_selection.reorg_en_calc(energies)
        print gen_result
        #matching up mol objects and energies
        mol_name_list = []
        name_list = cPickle.load( open(p_name, "rb" ) )
        for mol in name_list:
            temp_m = mol
            mol_name_list.append(temp_m)
        mol_en = []
        for i,mol in enumerate(mol_name_list):
           if mol.GetProp("_Energy") == "None":
               mol.SetProp("_Energy",gen_result[i][1])
               mol_en.append((mol,mol.GetProp("_Energy")))
        ml_sorted_by_en = sorted(mol_en, key=lambda tup: float(tup[1]))
        print ml_sorted_by_en
        #writing results pickle
        dump_list = []
        for mol,en in ml_sorted_by_en:
            name = mol.GetProp("_Name")
            en = mol.GetProp("_Energy")
            pm = AllChem.PropertyMol(mol)
            pm.SetProp("_Name", name)
            pm.SetProp("_Energy", en)
            dump_list.append(pm)
        cPickle.dump(dump_list, open(gen_name+"_results.p", "w+"))
        #setting crossover/elitism sizes
        crossover_size = int(args.population_size) * args.crossover_rate/100
        elite_size = int(args.population_size) * args.elitism_rate/100
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
        #print len(pairs_set)
        children = []
        for p in pairs_set:
            while True:
                try:
                    kids = ga_ring_crossover.ring_crossover(p[0],p[1])
                    children+=kids
                    break
                except (IndexError,RuntimeError,ValueError):
                    continue  
        #print len(children)
        #combining elite and crossed over mols
        mol_list_2 = []
        new_gen = elite_mols + children
        for i,nm in enumerate(new_gen):
            nm.SetProp("_Name", str(i)+"_test")
            mol_list_2.append(nm)
        #print len(mol_list_2)
        #print mol_list_2
        #removing duplicates and writing new generation
        new_mols = ga_setup_selection.check_previous_gens(mol_list_2)
        dump_list = []
        for m in new_mols:
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
        new_gen_name = str(args.population_size)+"_"+str(ng+1)+"_gen.p" 
        cPickle.dump(dump_list, open(new_gen_name, "w+"))
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
        com_files_neut = ga_setup_selection.com_generator_neutral_grndst(all_mol_coords)
        com_files_neg = ga_setup_selection.com_generator_neg_mol_neg_geom(all_mol_coords)
        com_files = com_files_neut + com_files_neg 
        job = group([run_gauss.s(com) for com in com_files])
        result = job.apply_async()
        result1 = result.get()
        opt_structs = [(name,struct) for name,en,struct in result1 if (name,struct)]
        first_energies = [(name,en) for name,en,struct in result1 if (name,en)]
        print first_energies
        neg_structs = [(name,struct) for name,struct in opt_structs if "neg" in name]
        neut_structs = [(name,struct) for name,struct in opt_structs if "neg" not in name]
        com_files_neg_2 = ga_setup_selection.com_generator_neg_mol_neut_charge(neg_structs)
        com_files_neut_2 = ga_setup_selection.com_generator_neut_mol_neg_charge(neut_structs)
        com_files_2 = com_files_neg_2 + com_files_neut_2
        job2 = group([run_gauss2.s(com) for com in com_files_2])
        result2 = job2.apply_async()
        second_energies = result2.get()
        print second_energies
        fe_sorted = sorted(first_energies, key=lambda tup: int(tup[0].split("_")[0]))
        se_sorted = sorted(second_energies, key=lambda tup: int(tup[0].split("_")[0]))
        energies = zip(fe_sorted,se_sorted)
        gen_result = ga_setup_selection.reorg_en_calc(energies)
        print gen_result
        #mol_name_list = []
        #name_list = cPickle.load( open(p_name, "rb" ) )
        #for mol in name_list:
        #    temp_m = mol
        #    mol_name_list.append(temp_m)
        mol_en = []
        i = 0
        for mol in mol_list:
           if mol.GetProp("_Energy") == "None":
               print mol.GetProp("_Name"), i
               mol.SetProp("_Energy",gen_result[i][1])
               i += 1
           mol_en.append((mol,mol.GetProp("_Energy")))
        ml_sorted_by_en = sorted(mol_en, key=lambda tup: float(tup[1]))
        print ml_sorted_by_en
        #writing results pickle
        dump_list = []
        for mol,en in ml_sorted_by_en:
            name = mol.GetProp("_Name")
            en = mol.GetProp("_Energy")
            pm = AllChem.PropertyMol(mol)
            pm.SetProp("_Name", name)
            pm.SetProp("_Energy", en)
            dump_list.append(pm)
        cPickle.dump(dump_list, open(gen_name+"_results.p", "w+"))
        #setting crossover/elitism sizes
        crossover_size = int(args.population_size) * args.crossover_rate/100
        elite_size = int(args.population_size) * args.elitism_rate/100
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
        #print len(pairs_set)
        children = []
        while True:
            try:
                kids = ga_ring_crossover.ring_crossover(p[0],p[1])
                children+=kids
                break
            except (IndexError,RuntimeError,ValueError):
                continue
        #print len(children)
        #combining elite and crossed over mols
        mol_list_2 = []
        new_gen = elite_mols + children
        for i,nm in enumerate(new_gen):
            nm.SetProp("_Name", str(i)+"_test")
            mol_list_2.append(nm)
        #print len(mol_list_2)
        #print mol_list_2
        #removing duplicates and writing new generation
        new_mols = ga_setup_selection.check_previous_gens(mol_list_2)
        dump_list = []
        for m in new_mols:
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
        new_gen_name = str(args.population_size)+"_"+str(ng+1)+"_gen.p"
        cPickle.dump(dump_list, open(new_gen_name, "w+"))

