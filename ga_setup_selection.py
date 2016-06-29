from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import random
import cPickle
import glob

class MB():
    '''simple class that is used for checking molecules are the same using the daylight fingerprint and Tanimoto similarity'''

    def __init__(self, mol):
        self.mol = mol
        self.bitstring = FingerprintMols.FingerprintMol(mol)
    
    def __eq__(self,other):
        return DataStructs.FingerprintSimilarity(self.bitstring,other.bitstring) == 1.0




def chunks(l, n):
    '''makes a list into n sized chunks'''
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
              "%nprocshared=6",
              "%chk=molecule_0",
              "#B3LYP/6-311+G** FChk NoSymmetry Opt=(ModRedundant,)",
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
              "%nprocshared=6",
              "%chk=molecule_0",
              "#B3LYP/6-311+G** FChk NoSymmetry Opt=(ModRedundant,)",
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
               "%nprocshared=6",
               "%chk=molecule_0",
               "#B3LYP/6-311+G** FChk NoSymmetry",
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
              "%nprocshared=6",
              "%chk=molecule_0",
              "#B3LYP/6-311+G** FChk NoSymmetry",
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
    results_ps = glob.glob("*_results.p")
    prev_mols_list = []
    for re in results_ps:
        ml_list = cPickle.load( open(re, "rb" ) )
        for ml in ml_list:
            prev_mols_list.append(ml)
    print "previous molecules to check "+str(len(prev_mols_list))
    smiles_en_list = [(Chem.MolToSmiles(ml),(ml.GetProp("_Fit"),ml.GetProp("_Energy"),ml.GetProp("_Ea"))) for ml in prev_mols_list] 
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
                    mol.SetProp("_Energy", energy[1])
                    mol.SetProp("_Fit", energy[0])
                    mol.SetProp("_Ea", energy[2])
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

def check_previous_gens_fp(new_mol_list):
    '''similar function to above but using molecule fps instead of smiles strings
    should lead to fewer resonance structures making it through'''
    results_ps = glob.glob("*_results.p")
    print results_ps
    prev_mols_list = []
    for re in results_ps:
        ml_list = cPickle.load( open(re, "rb" ) )
        for ml in ml_list:
            prev_mols_list.append(ml)
    print "previous mols to check "+str(len(prev_mols_list))
    uniq_prev_mols = list()
    for pm in prev_mols_list:
        pmb = MB(pm)
        if pmb not in uniq_prev_mols:
            uniq_prev_mols.append(pmb)
    print "unique mols "+str(len(uniq_prev_mols))
    matches = []
    for mol in new_mol_list:
        nmb = MB(mol)
        for pmb in uniq_prev_mols:
            if nmb == pmb:
                pm = pmb.mol
                fit = pm.GetProp("_Fit")
                en = pm.GetProp("_Energy")
                ea = pm.GetProp("_Ea")
                mol.SetProp("_Fit", fit)
                mol.SetProp("_Energy", en)
                mol.SetProp("_Ea", ea)
                matches.append(mol)
    print "number of matches "+str(len(matches))
    new_matches = matches + new_mol_list
    new_matches_set = set(new_matches)
    new_matches_list = list(new_matches_set)
    uniq_ng = len(new_mol_list) - len(matches)
    print "unqiue mols in gen "+str(uniq_ng)
    new_matches_sorted = sorted(new_matches_list, key=lambda x: (int(x.GetProp("_Name").split("_")[0]))) 
    return new_matches_sorted


def init_gen_dumper(p_name,mol_list):
    '''dumps clean generation to pickle'''
    dump_list = []
    for m in mol_list:
        name = m.GetProp("_Name")
        pm = AllChem.PropertyMol(m)
        pm.SetProp("_Name", name)
        pm.SetProp("_Energy", None)
        pm.SetProp("_Ea", None)
        pm.SetProp("_Fit", None)
        dump_list.append(pm)
    cPickle.dump(dump_list, open(p_name, "w+"))

def gen_loader(p_name):
    '''loads up a generation from pickle'''
    mol_name_list = []
    name_list = cPickle.load( open(p_name, "rb" ) )
    for mol in name_list:
        temp_m = mol
        mol_name_list.append(temp_m)
    return mol_name_list

def results_dumper(p_name,ml_sorted_by_en):
    '''dumps gen results to pickle'''
    gen_name = p_name.split(".")[0]
    dump_list = []
    for data in ml_sorted_by_en:
        mol = data[0]
        name = mol.GetProp("_Name")
        en = mol.GetProp("_Energy")
        ea = mol.GetProp("_Ea")
        fit = mol.GetProp("_Fit")
        pm = AllChem.PropertyMol(mol)
        pm.SetProp("_Name", name)
        pm.SetProp("_Fit", fit)
        pm.SetProp("_Energy", en)
        pm.SetProp("_Ea", ea)
        dump_list.append(pm)
    cPickle.dump(dump_list, open(gen_name+"_results.p", "w+"))

def new_gen_dumper(pop_size,ng,new_mols):
    '''dumps new gen to pickle, keeping energies from matched molecules'''
    dump_list = []
    for m in new_mols:
        name = m.GetProp("_Name")
        try:
            energy = m.GetProp("_Energy")
        except KeyError:
            en = "None"
            energy = m.SetProp("_Energy", en)
        try:
            ea = m.GetProp("_Ea")
        except KeyError:
            eleca = "None"
            ea = m.SetProp("_Ea", eleca)
        try:
            fit = m.GetProp("_Fit")
        except KeyError:
            fitness = "None"
            fit = m.SetProp("_Fit", fitness)
        pm = AllChem.PropertyMol(m)
        pm.SetProp("_Name", name)
        pm.SetProp("_Energy", energy)
        pm.SetProp("_Fit", fit)
        pm.SetProp("_Ea", ea)
        dump_list.append(pm)
    new_gen_name = str(pop_size)+"_"+str(ng+1)+"_gen.p"
    cPickle.dump(dump_list, open(new_gen_name, "w+"))

def result_matcher(gen_result,mol_name_list,maximise):
    '''matches results from calcs with correct mols'''
    mol_en = []
    print gen_result
    print mol_name_list
    i = 0
    for mol in mol_name_list:
        if mol.GetProp("_Energy") == "None":
            try:
                print mol.GetProp("_Name"),gen_result[i]
                mol.SetProp("_Energy",gen_result[i][1])
                mol.SetProp("_Ea", gen_result[i][2])
                i += 1
            except IndexError:
                print "results list wrong length"
        mol_en.append((mol,mol.GetProp("_Energy"), mol.GetProp("_Ea")))
    print mol_en
    mol_fit = []
    k1 = 1.0
    k2 = 1.0
    for mol in mol_en:
        reorg = float(mol[1])
        ea = float(mol[2])
        if ea < 3.0:
           k2 = 3 - ea
        elif ea > 4.0:
           k2 = ea - 4
        elif 3 < ea < 4:
           k2 = 0
        fit = k1*reorg + k2*ea
        mol[0].SetProp("_Fit", fit)
        mol_fit.append((mol[0], fit, ea, reorg))
    print maximise
    if maximise == "True":
        print " max = True"
        ml_sorted_by_en = sorted(mol_en, key=lambda tup: float(tup[1]),reverse=True)
    if maximise == "False":
        print "max = False"
        ml_sorted_by_en = sorted(mol_fit, key=lambda tup: float(tup[1]))
    return ml_sorted_by_en

def make_pairs(crossover_mols):
    '''makes pairs from ranked and selected mols, pairs first to last and so on'''
    pairs_set = []
    pairs = zip(crossover_mols[::1], crossover_mols[-1::-1])
    for p in pairs:
        if (p[0],p[1]) in pairs_set:
            pass
        if (p[1],p[0]) in pairs_set:
            pass
        else:
            pairs_set.append((p[0],p[1]))
    return pairs_set

def write_overall_results(pop_size):
    '''writes one pickle containing all unique molecules sampled during the search'''
    overall = []
    result_ps = glob.glob("*_results.p")
    for rp in result_ps:
        gen = rp.split("_")[1]
        temp_list = []
        mr_list = cPickle.load( open(rp, "rb" ) )
        for mol in mr_list:
            temp_list.append((mol,gen))
        overall+=temp_list
    print "total molecules sampled "+str(len(overall))
    overall_rang = len(overall)
    overall_sorted = sorted(overall, key=lambda tup: float(tup[0].GetProp("_Fit")))
    for o in overall_sorted:
        mol = o[0]
        gen = o[1]
        mol.SetProp("_gen", gen)
        #bs = FingerprintMols.FingerprintMol(mol)
        mol2 = MB(mol)
        #print mol2.bitstring
        #mol_set.add(mol2)
        if mol2 not in mol_set:
            mol_set.append(mol2)
    print "total unique molecules "+str(len(mol_set))
    final_mol_list = []
    for m in mol_set:
        mol = m.mol
        gen = mol.GetProp("_gen")
        name = mol.GetProp("_Name")
        new_name = name+"_"+gen
        mol.SetProp("_Name", new_name)
        final_mol_list.append(mol)
    cPickle.dump(final_mol_list, open(str(pop_size)+"_final.p", "w+"))

