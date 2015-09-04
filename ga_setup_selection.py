from rdkit import Chem
from rdkit.Chem import AllChem
import random

def chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]

def mutator(smiles_list, mutation_factor, N_number): 
    '''works on the list form of a SMILES string, loops over every atom, with mutation_factor chance to mutate, 
    loop continues until N_number of nitrogen atoms are present'''
    ia = 0 #counter for getting right index in bond_list
    #print smiles_list
    while smiles_list.count('n') != N_number: 
        smiles_string = ''.join(smiles_list) #need to change smiles_list back into string 
        m_test = Chem.MolFromSmiles(smiles_string) #to change back into a molecule
        bond_list = [len(a.GetBonds()) for a in m_test.GetAtoms()] #to get number of bonds for each atom
        for i,atom in enumerate(smiles_list):
            try:
                int(atom) #integers are present in the smiles string so these need to be ignored
            except ValueError: #only items which fail to turn into an integer will be mutated
                try:
                    if int(smiles_list[i+1]) == True:
                        pass
                except ValueError:
                    if bond_list[ia] > 2: #skip central atoms
                        pass
                    elif random.randint(0,100) <= mutation_factor: #if randint is within mutation_factor 
                        if smiles_list.count('n') < N_number:#and not enough N
                            atom = 'n' #mutate
                        elif smiles_list.count('n') > N_number: #if too many N
                            atom = 'c' #mutate back to C
                    smiles_list[i]=atom #put back in the list correctly
                    ia += 1
                    if ia >= len(bond_list): #stops loop breaking 
                        ia = 0
        test_string = ''.join(smiles_list)
        test_m = Chem.MolFromSmiles(test_string)
        if test_m == None:
            test_list = list(test_string)
            mutator(test_list,mutation_factor = 50, N_number = N_number)
        if smiles_list.count('n') == N_number: 
            return smiles_list #returns list of atoms

def generate_initial_population(smiles_string, pop, N_number):
    '''generates an intial mutated population based on one molecule, including N_number nitrogens 
    in each member upto a pop sized set'''
    population = set() #sets are unordered collections of unique objects, should stop identical ones from slipping in
    while len(population) < pop: #build up to pop size
        smiles_list = list(smiles_string) #converting to list for mutator
        new_member_list = mutator(smiles_list, mutation_factor = 10, N_number = N_number) 
        new_member = ''.join(new_member_list)
        if new_member not in population: #if string isn't in set
            if new_member[::-1] not in population: #and reverse of the string 
                population.add(new_member) #add it, rdkit has mol comparison tools for when this needs sophistication
    if len(population) == pop: #end when pop size reached
        return population
    
def mol_name_maker(pop):
    '''adds names on to the pop list as a property of the molecule, mol_list is easier to work with need to be careful
    adding to it however'''
    i = 1 #just a counter for the name
    mol_list = []
    for p in pop:
        temp_m = Chem.MolFromSmiles(p) #builds mol from smiles
        temp_m.SetProp("_Name",str(i)+"_test") #setting some generic name
        mol_list.append(temp_m) #add to the list
        i += 1
    return mol_list

def coords_generator(mol_list):
    '''generates a set of 3D coordinates for each molecule in mol_list, returns the entire list'''
    all_mol_coords = []
    for m in mol_list:
        #pm = Chem.MolFromSmiles(p)
        #print m.GetProp("_Name")  
        if m is None: #test to make sure molecule is the correct object
            print "Molecule is junk"
        pmH = Chem.AddHs(m) #mols are represented sans H, this adds them back in
        AllChem.EmbedMolecule(pmH) #needs to be embedded for 3D coords
        AllChem.UFFOptimizeMolecule(pmH) #cheap optimisation step, improves conformation and later energies
        coords = Chem.MolToMolBlock(pmH) #this is a mol block of the coords after opt
        coords_list = []
        for line in coords.split("\n"):  #looping over each line ignoring everything except atom labels and coords
            test =line.strip().split()
            if len(test) == 0:
                pass
            else:
                try:
                    if int(test[0]):
                        pass
                except ValueError:
                    if len(test) > 2:
                        label = test[3]
                #tmp.append(label)
                        coords_list.append((label,test[0:3])) 
        all_mol_coords.append(coords_list) #add each mols coords to overall coords list
    return all_mol_coords
    
    
def com_generator_neutral_grndst(coords_list):
    ''''write comfile for neutral ground state calc for a list of coordinates, finishes when all are written'''
    all_mol_files = []
    i=1
    for mol in coords_list:
        mol_file = []
        name = str(i)+"_test" #set name, should still be in same order as original mol_list
        #com_f = open(name+".com", "w")
        gs = ["%mem=1GB",
              "%nprocshared=8",
              "%chk=molecule_0",
              "#B3LYP/6-31G** FChk NoSymmetry Opt=(ModRedundant,)",
              "", 
              "neutral minimisation "+name,
              "",
              "0 1"]
        for line in gs:
             mol_file.append(line+"\n")
        #     com_f.write(line+"\n")
        for atom in mol:
            if "_test" in atom:
                pass
            else:
                mol_file.append(atom)
         #       com_f.write("%6s %13.6f %13.6f %13.6f\n"% (atom[0],float(atom[1][0]),float(atom[1][1]),float(atom[1][2])))
        #com_f.write(3*"\n")
        #com_f.close()
        i +=1
        all_mol_files.append(mol_file)
    return all_mol_files
    print "Com files written"
    
def com_generator_neg_mol_neg_geom(coords_list):
    '''write com file for negative mol geom opt'''
    all_mol_files = []
    i = 1
    for mol in coords_list:
        mol_file = []
        #com_f = open(str(i)+"_test_neg.com", "w")
        gs = ["%mem=1GB",
              "%nprocshared=8",
              "%chk=molecule_0",
              "#B3LYP/6-31G** FChk NoSymmetry Opt=(ModRedundant,)",
              "", 
              "negative minimisation "+str(i)+"_test_neg",
              "",
              "-1 2"]
        for line in gs:
             mol_file.append(line+"\n")
             #com_f.write(line+"\n")
        for atom in mol:
            #com_f.write("%6s %13.6f %13.6f %13.6f\n"% (atom[0],float(atom[1][0]),float(atom[1][1]),float(atom[1][2])))
            mol_file.append(atom)
        #com_f.write(3*"\n")
        #com_f.close()
        i += 1
        all_mol_files.append(mol_file)
    return all_mol_files
    print "Negative mol, negative geom com files written"

def com_generator_neg_mol_neut_charge(coords_list):
    '''write com files for opt structures with changed multiplicities'''
    all_mol_files = []
    for mol in coords_list:
        struct = mol[1]
        name = mol[0]+"_opt"
        mol_file = []
        gs = ["%mem=1GB",
               "%nprocshared=8",
               "%chk=molecule_0",
               "#B3LYP/6-31G** FChk NoSymmetry",
               "",
               "single point energy "+name,
               "",
               "0 1"]
        for line in gs:
            mol_file.append(line+"\n")
        for atom in struct:
            mol_file.append(atom)
        all_mol_files.append(mol_file)
    return all_mol_files

def com_generator_neut_mol_neg_charge(coords_list):
    '''write com files for opt structures with changed multiplicities'''
    all_mol_files = []
    for mol in coords_list:
        struct = mol[1]
        name = mol[0]+"_opt"
        mol_file = []
        gs = ["%mem=1GB",
              "%nprocshared=8",
              "%chk=molecule_0",
              "#B3LYP/6-31G** FChk NoSymmetry",
              "",
              "single point energy "+name,
              "",
              "-1 2"]
        for line in gs:
            mol_file.append(line+"\n")
        for atom in struct:
            mol_file.append(atom)
        all_mol_files.append(mol_file)
    return all_mol_files


def roulette_wheel_selection(mol_list_en,elitism_rate,crossover_rate):
    '''roulette wheel selection of mols to be included in elitism and crossover'''
    elite_mols = []
    crossover_mols = []
    total_fitness = sum(float(mol[1]) for mol in mol_list_en)
    norm_mol_fitness = []
    for mol,en in mol_list_en:
        norm_fitness = float(en)/total_fitness
        norm_mol_fitness.append((mol,norm_fitness))
    fittest_mols = []
    crossover_mols = []
    norm_max = sum(float(mol[1]) for mol in norm_mol_fitness)
    rel_fitness = [norm_fit for mol,norm_fit in norm_mol_fitness if norm_fit]
    wrong_probs = [sum(rel_fitness[:i+1]) for i in range(len(rel_fitness))] 
    probs = [new_max - wp for wp in wrong_probs]
    for i,mol in enumerate(norm_mol_fitness):
        pick = random.uniform(0, new_max) #random number to decide what happens to mol
        if probs[i] > pick: #if a mols fitness is higher than the rand number add it to untouched
            if len(fittest_mols) >= elitism_rate: #only pass a certain % through
                pass
            else:
                fittest_mols.append(mol)
                print "pocket size = "+str(probs[i])
                print "rand number = "+str(pick)
    print len(fittest_mols)
    for i,mol in enumerate(mol_fitness):
        while len(crossover_mols) < 200:
            pick = random.uniform(0,new_max)
            if probs[i] > pick:
                crossover_mols.append(mol)
                print "pocket_size = "+str(probs[i])
                print "rand number = "+str(pick)
    print len(crossover_mols)
    return fittest_mols, crossover_mols

def sus_selection(mol_list_en, elitism_rate, new_pop_size):
    '''related to roulette wheel selection except many equally sized steps are taken around the wheel
    elitism is included by taking the first elitism_rate members of the entire new selected pop'''
    from numpy import arange
    elite_mols = []
    crossover_mols = []
    total_fitness = sum(float(mol[1]) for mol in mol_list_en)
    norm_fitness = [float(en)/total_fitness for mol,en in mol_list_en]
    norm_total = sum(nf for nf in norm_fitness)
    wrong_probs = [sum(norm_fitness[:i+1]) for i in range(len(norm_fitness))] 
    probs = [norm_total - wp for wp in wrong_probs]
    mol_only = [mol for mol,en in mol_list_en if mol]
    mol_probs = zip(mol_only,probs)
    step_size = norm_total/new_pop_size
    #print step_size
    start_point = random.uniform(0,step_size)
    print start_point
    points = arange(start_point,(start_point+(step_size*new_pop_size)),step_size)
    print len(points)
    if len(points) > new_pop_size:
       points = arange(start_point,(start_point+(step_size*new_pop_size)),step_size) 
    for i,step in enumerate(points):
        try:
            if mol_probs[i][1] > step:
                crossover_mols.append(mol_probs[i][0])
        except IndexError:
            pass
    for i,mol in enumerate(crossover_mols):
        if i < elitism_rate:
            elite_mols.append(mol)
    return crossover_mols, elite_mols
    
def tournament_selection(mol_list_en, elitism_rate, new_pop_size):
    '''tournament selection pits two randomly chosen members of the population against one another,
    the fittest has a probability p to win, the weaker 1-p'''
    elite_mols = []
    crossover_mols = []
    p0 = 75
    p1 = 76
    for i,mol in enumerate(mol_list_en):
        if i < elitism_rate:
            elite_mols.append(mol[0])
    total_fitness = sum(float(mol[1]) for mol in mol_list_en)
    proper_fitness = [total_fitness - float(mol[1]) for mol in mol_list_en]
    while len(crossover_mols) < new_pop_size:
        challenger_1 = random.randint(0,(new_pop_size-1))
        challenger_2 = random.randint(0,(new_pop_size-1))
        print str(challenger_1)+" versus "+str(challenger_2)
        c1_fit = proper_fitness[challenger_1]
        c2_fit = proper_fitness[challenger_2]
        winner = random.randint(0,100)
        if winner <= p0:
            print "stronger won"
            crossover_mols.append(mol_list_en[challenger_1][0])
        if winner > p1:
            print "weaker won"
            crossover_mols.append(mol_list_en[challenger_2][0])
    return crossover_mols,elite_mols
    
def rank_selection(mol_list_en, elitism_rate, new_pop_size):
    '''similar to roulette wheel selection but relative fitness is proportional to rank instead of
    absolute fitness, should help improve diversity and smooth out large fitness differences'''
    elite_mols = []
    crossover_mols = []
    mol_sorted_by_en = sorted(mol_en, key=lambda tup: float(tup[1]))
    mol_rank = [(i,mol[0]) for i, mol in enumerate(mol_sorted_by_en)]
    elite_mols = [mol for i,mol in mol_rank if i < 50]
    total_rank = sum(i for i,mol in mol_rank)
    norm_rank = [float(rank)/total_rank for rank,mol in mol_rank]
    norm_total = sum(nf for nf in norm_rank)
    wrong_probs = [sum(norm_rank[:i+1]) for i in range(len(norm_rank))] 
    probs = [norm_total - wp for wp in wrong_probs]
    while len(crossover_mols) < new_pop_size:
        spin = random.randint(0,(len(mol_rank)-1))
        pick = random.uniform(0, norm_total)
        if probs[spin] > pick:
            crossover_mols.append(mol_rank[spin][1])
    return crossover_mols,elite_mols      

def check_previous_gens(new_mol_list):
    '''checks previous results.p for identical molecules from crossover'''
    import glob
    import cPickle
    unique_mols = []
    results_ps = glob.glob("*_results.p")
    prev_mols_list = []
    for re in results_ps:
        ml_list = cPickle.load( open(re, "rb" ) )
        for ml in ml_list:
            prev_mols_list.append(ml)
    print "previous molecules to check "+str(len(prev_mols_list))
    smiles_en_list = [(Chem.MolToSmiles(ml),ml.GetProp("_Energy")) for ml in prev_mols_list]
    smiles_en_dict = dict(smiles_en_list)
    unique_se = [(k,v) for k,v in smiles_en_dict.iteritems()]
    print "number of previous unique molecules "+str(len(unique_se))
    #print test_3
    matches = []
    for mol in new_mol_list:
        i = 0
        for msmiles, energy in unique_se:
            if Chem.MolToSmiles(mol) == msmiles:
                i += 1
                if i > 1:
                    break
                else:
                    mol.SetProp("_Energy", energy)
                    matches.append(mol)
            if i > 1:
                break
    print "number of matches "+str(len(matches))
    new_matches = matches + new_mol_list
    new_matches_set = set(new_matches)
    uniq_ng = len(new_mol_list) - len(matches)
    print "unqiue mols in gen "+str(uniq_ng)
    new_matches_list = list(new_matches_set)
    new_matches_sorted = sorted(new_matches_list, key=lambda x: (int(x.GetProp("_Name").split("_")[0])))
    return new_matches_sorted 

def reorg_en_calc(energies):
    '''calculate reorganisation energies from run results'''
    result = []
    sep_lists = chunks(energies, 2)
    for molen in sep_lists:
        print molen
        x = (float(molen[1][1][1]) - float(molen[1][0][1])) 
        y =  (float(molen[0][1][1]) - float(molen[0][0][1]))
        re_ha = x + y
        re_ev = re_ha *27.212
        result.append((molen[0][0][0], re_ev))
    return result
