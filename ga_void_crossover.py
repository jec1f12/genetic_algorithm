import random
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

##crossover for void GA that will probably not make it 


class MB():#class uses python eq method, if the fingerprint == 1 the molecules are identical
    '''simple class that is used for checking molecules are the same using the daylight fingerprint and Tanimoto similarity'''

    def __init__(self, mol):
        self.mol = mol
        self.bitstring = FingerprintMols.FingerprintMol(mol)

    def __eq__(self,other):
        return DataStructs.FingerprintSimilarity(self.bitstring,other.bitstring) == 1.0



def mutator_void(mol):
    '''function for mutating a biphenyl through the addition of substituents outside the ring system, returns a mutated mol'''
    import random
    atom_sites = [0,4,9,11,14,17,18,19,20,21]#these are the atom numbers for the ring substituents
    subs_list = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,7,7,7,6,6,6,6,6,6,6,9,17,35,53]#atomic number of substituents
    carbon_subs_smiles = ['C[Si](C)(C)(C)','C(C)(C)(C)','C[Si](CC)(CC)(CC)','C(F)(F)(F)'
                      ,'C#N','C#C', 'C']#smiles for carbon subs
    nitrogen_subs_smiles = ['*[H]','*(C)(C)','*(=O)(=O)']#nitrogen subs
    pmH = Chem.AddHs(mol)#this ensures hydrogens are there so all sites can be mutated
    for atom in atom_sites:
        replacement = random.choice(subs_list)
        print atom, replacement
        if replacement == 6:#for carbon and nitrogen specific bonding must be used for various subs
            print "carbon"
            carbon_sub = random.choice(carbon_subs_smiles)
            if carbon_sub == 'C':#methyl is easiest
                pmHe = Chem.EditableMol(pmH)
                pmHe.ReplaceAtom(atom, Chem.Atom(6))
                tm2 = pmHe.GetMol()
                Chem.SanitizeMol(tm2)
                pmH = tm2
            else:
                temp_mol = Chem.MolFromSmiles(carbon_sub)
                comb_mol = Chem.CombineMols(pmH,temp_mol)
                frags = Chem.GetMolFrags(comb_mol)
                pmHe = Chem.EditableMol(comb_mol)
                pmHe.ReplaceAtom(atom, Chem.Atom(6))
                if '[Si]' in carbon_sub:
                    pmHe.AddBond(atom,frags[1][0],Chem.BondType.TRIPLE)#bond from first atom in subs to atomic site
                    tm2 = pmHe.GetMol()
                    Chem.SanitizeMol(tm2)
                    pmH = tm2
                elif '#' in carbon_sub:
                    pmHe.AddBond(atom,frags[1][1],Chem.BondType.TRIPLE)#bond from atom to second atom in subs
                    pmHe.RemoveAtom(frags[1][0])
                    tm2 = pmHe.GetMol()
                    Chem.SanitizeMol(tm2)
                    pmH = tm2
                else:
                    pmHe.AddBond(atom,frags[1][1],Chem.BondType.SINGLE)#all bonds need to be specified for tetra subs
                    pmHe.AddBond(atom,frags[1][2],Chem.BondType.SINGLE)
                    pmHe.AddBond(atom,frags[1][3],Chem.BondType.SINGLE)
                    pmHe.RemoveAtom(frags[1][0])
                    tm2 = pmHe.GetMol()
                    Chem.SanitizeMol(tm2)
                    pmH = tm2
        if replacement == 7:
            print "nitrogen"
            nitro_sub = random.choice(nitrogen_subs_smiles)
            print nitro_sub
            temp_mol = Chem.MolFromSmiles(nitro_sub)
            comb_mol = Chem.CombineMols(pmH,temp_mol)
            frags = Chem.GetMolFrags(comb_mol)
            pmHe = Chem.EditableMol(comb_mol)
            pmHe.ReplaceAtom(atom,Chem.Atom(7))
            if "O" in nitro_sub:#like carbon bonding rules need to be followed
                pmHe.AddBond(atom,frags[1][1],Chem.BondType.DOUBLE)#specifically rdkits NO2 bonding info
                pmHe.AddBond(atom,frags[1][2],Chem.BondType.DOUBLE)#this does result in the correct ionic form once sanitised
                pmHe.RemoveAtom(frags[1][0])
                tm2 = pmHe.GetMol()
                Chem.SanitizeMol(tm2)
                pmH = tm2
            elif "C" in nitro_sub:
                pmHe.RemoveAtom(frags[1][0])
                pmHe.AddBond(atom,frags[1][0],Chem.BondType.SINGLE)
                pmHe.AddBond(atom,frags[1][1],Chem.BondType.SINGLE)
                tm2 = pmHe.GetMol()
                Chem.SanitizeMol(tm2)
                pmH = tm2
            else:
                pmHe.ReplaceAtom(atom,Chem.Atom(7))
                pmHe.RemoveAtom(frags[1][0])
                tm2 = pmHe.GetMol()
                Chem.SanitizeMol(tm2)
                pmH = tm2
        else:
            pmHe = Chem.EditableMol(pmH)
            pmHe.ReplaceAtom(atom,Chem.Atom(replacement))#simplest for single atom subs just replace
            tm2 = pmHe.GetMol()
            Chem.SanitizeMol(tm2)
            pmH = tm2
    return tm2
    
    
def generate_initial_population_void(mol,pop_size):
    '''generates an initial biphenyl population and checks for uniqueness using molecular fps'''
    population_init = []
    population = []
    while len(population_init) < pop_size:
        new_member = mutator_void(mol)
        new_mem_fp = MB(new_member)
        if new_mem_fp not in population_init:
            population_init.append(new_mem_fp)
    for p in population_init:
        pmol = p.mol
        population.append(pmol)
    return population

def get_subs(mol):
    '''gets subs list from a biphenyl, list returned contains atomic site, atom type and mol object of sub'''
    molh = Chem.AddHs(mol)
    rings = molh.GetRingInfo()
    ring_ats = rings.AtomRings()
    flat_ring_ats =  sum(ring_ats, ())#just makes a flat list out of the pair of tuples returned by ring.AtomRings()
    list_atom_sites = []
    for ring_at in flat_ring_ats:
        n = [n.GetIdx() for n in molh.GetAtomWithIdx(ring_at).GetNeighbors()]#get neighbours for all our ring atoms
        non_ring_n = [nrn for nrn in n if nrn not in flat_ring_ats]#gets non ring neighbours
        list_atom_sites.append(non_ring_n)
    actual_atom_sites = [atom for sublist in list_atom_sites for atom in sublist]#gets a list of all our substituent atoms
    print actual_atom_sites
    subs_list = []
    for at_s in actual_atom_sites:
        for atoms in molh.GetAtoms():
            if atoms.GetIdx() == at_s:
                #print atoms.GetSymbol()
                n = [n.GetIdx() for n in atoms.GetNeighbors()]
                #print n
                mol_e = Chem.EditableMol(mol)
                mol_e.RemoveBond(n[0],atoms.GetIdx())
                mol_e.ReplaceAtom(atoms.GetIdx(),Chem.Atom(0))#dummy atoms allow us to pass sanitisation steps 
                broke_mol = mol_e.GetMol()
                Chem.SanitizeMol(broke_mol)
                frags = Chem.GetMolFrags(broke_mol,asMols=True)
                if at_s == 0:
                    subs_list.append((atoms.GetIdx(),atoms.GetAtomicNum(),frags[0]))
                else:
                    subs_list.append((atoms.GetIdx(),atoms.GetAtomicNum(),frags[1]))
    return subs_list

def non_ring_atoms(mol):
    '''function that returns the non ring atoms of a molecule'''
    rings = mol.GetRingInfo()
    ring_ats = rings.AtomRings()
    test_2 =  sum(ring_ats, ())
    #print test_2
    atomz = [atom.GetIdx() for atom in mol.GetAtoms()]
    #not_ir_ats = [item for sublist in actual_cp for item in sublist[1] if item not in ring_ats]
    #x for x in t if x not in s
    not_ir_ats = [at for at in atomz if at not in test_2]
    #print not_ir_ats
    list_atom_sites = []
    for ring_at in test_2:
    #print ring_at
        n = [n.GetIdx() for n in mol.GetAtomWithIdx(ring_at).GetNeighbors()]
    #print n
        non_ring_n = [nrn for nrn in n if nrn not in test_2]
        list_atom_sites.append(non_ring_n)
    actual_atom_sites = [atom for sublist in list_atom_sites for atom in sublist]
    #print actual_atom_sites
    sorted_aas = sorted(actual_atom_sites)
    return sorted_aas

def add_subs2(mol,subs_list):
    '''function that adds all the subs in subs_list to a molecule'''
    i = 0
    blank_molh = Chem.AddHs(mol)
    sort_aas = non_ring_atoms(blank_molh)
    print sort_aas
    for aas in sort_aas:
        #print "old non ring atom "+str(aas)
        subs = subs_list[i]
        #print subs
        sub_mol = subs[2]
        sub_mol_smi = Chem.MolToSmiles(sub_mol)
        comb_mol = Chem.CombineMols(blank_molh,sub_mol)
        frags = Chem.GetMolFrags(comb_mol)
        #print frags
        comb_mole = Chem.EditableMol(comb_mol)
        n = [n.GetIdx() for n in comb_mol.GetAtomWithIdx(aas).GetNeighbors()]
        print n
        if "O=" in sub_mol_smi:
            comb_mole.AddBond(frags[1][0],n[0],Chem.BondType.SINGLE)
            comb_mole.RemoveBond(frags[1][0],frags[1][2])
            comb_mole.RemoveBond(frags[1][0],frags[1][1])
            comb_mole.AddBond(frags[1][0],frags[1][1],Chem.BondType.DOUBLE)
            comb_mole.AddBond(frags[1][0],frags[1][2],Chem.BondType.DOUBLE)
            comb_mole.ReplaceAtom(frags[1][0],Chem.Atom(subs[1]))
            comb_mole.ReplaceAtom(aas, Chem.Atom(0))
            comb_mole.ReplaceAtom(n[0], Chem.Atom(0))
            cm2 = comb_mole.GetMol()
            Chem.SanitizeMol(cm2)
            blank_molh = cm2 
            i += 1  
        else:
            comb_mole.AddBond(frags[1][0],n[0],Chem.BondType.SINGLE)
            comb_mole.ReplaceAtom(frags[1][0],Chem.Atom(subs[1]))
            comb_mole.ReplaceAtom(aas, Chem.Atom(0))
            comb_mole.ReplaceAtom(n[0], Chem.Atom(0))
            cm2 = comb_mole.GetMol()
            Chem.SanitizeMol(cm2)
            blank_molh = cm2 
            i += 1 
    return blank_molh

def dummy_atom_remover_non_ring(mol):
    '''removing dummy atoms added in by subs adder is a two step process, non rings are replaced afterwards'''
    nring_dum_ats = []
    for atoms in mol.GetAtoms():
        if atoms.GetAtomicNum() == 0:
            n = [n.GetIdx() for n in atoms.GetNeighbors()]
            if len(n) < 2:
                #print "non ring dummy atom "+str(atoms.GetIdx())
                nring_dum_ats.append(atoms.GetIdx())
    #print nring_dum_ats
    snring_dum_ats = sorted(nring_dum_ats,reverse=True)
    print snring_dum_ats
    for nrats in snring_dum_ats:
        mol_e = Chem.EditableMol(mol)
        mol_e.RemoveAtom(nrats)
        mol2 = mol_e.GetMol()
        Chem.SanitizeMol(mol2)
        mol = mol2
    for atoms in mol.GetAtoms():
        if atoms.GetAtomicNum() == 0:
            print atoms.GetIdx()
            mol_e = Chem.EditableMol(mol)
            mol_e.ReplaceAtom(atoms.GetIdx(),Chem.Atom(6))
            mol2_nonsan = mol_e.GetMol()
            mol = mol2_nonsan
    return mol

def dummy_atom_remover_ring(mol):
    '''this replaces all the dummy atoms inside the biphenyl rings'''
    rings = mol.GetRingInfo()
    ring_ats = rings.AtomRings()
    flat_rats =  sum(ring_ats, ())
    #print test_2
    #print ring_dum_ats
    for atoms in mol.GetAtoms():
        if atoms.GetIdx() in flat_rats:
            print atoms.GetIdx(), atoms.GetIsAromatic()
            if atoms.GetIsAromatic() == False:
                atoms.SetIsAromatic(True)
    Chem.SanitizeMol(mol)
    return mol
    
def void_single_point(mol1,mol2):
    import random
    children = []
    n = random.randint(0,9)
    print n
    mol1_subs_list = get_subs(mol1)
    mol2_subs_list = get_subs(mol2)
    print mol1_subs_list
    print mol2_subs_list
    child1_subs_list = mol1_subs_list[:n] + mol2_subs_list[n:]
    child2_subs_list = mol2_subs_list[:n] + mol1_subs_list[n:]
    blank_mol = Chem.MolFromSmiles('Clc1cc(Cl)c(-c2cc(Cl)c(Cl)cc2Cl)cc1Br')
    child_mol1_dum = add_subs2(blank_mol, child2_subs_list)
    child_mol2_dum = add_subs2(blank_mol, child1_subs_list)
    child_mol1_nring_dum = dummy_atom_remover_non_ring(child_mol1_dum)
    child_mol2_nring_dum = dummy_atom_remover_non_ring(child_mol2_dum)
    child_mol1 = dummy_atom_remover_ring(child_mol1_nring_dum)
    child_mol2 = dummy_atom_remover_ring(child_mol2_nring_dum)
    children.append(child_mol1)
    children.append(child_mol2)
    return children

def void_double_point(mol1,mol2):
    import random
    children = []
    n1 = random.randint(0,8)
    print n1
    n2 = random.randint((n1+1),9)
    print n2
    mol1_subs_list = get_subs(mol1)
    mol2_subs_list = get_subs(mol2)
    print len(mol1_subs_list)
    print len(mol2_subs_list)
    child_list_1 = mol1_subs_list[:n1] + mol2_subs_list[n1:n2]
    child_list_2 = mol2_subs_list[:n1] + mol1_subs_list[n1:n2]
    print child_list_1
    print child_list_2
    child_full_list1 = child_list_1 + mol1_subs_list[len(child_list_1):]
    child_full_list2 = child_list_2 + mol2_subs_list[len(child_list_2):]
    print child_full_list1, len(child_full_list1)
    print child_full_list2, len(child_full_list2)
    blank_mol = Chem.MolFromSmiles('Clc1cc(Cl)c(-c2cc(Cl)c(Cl)cc2Cl)cc1Br')
    child_mol1_dum = add_subs2(blank_mol, child_full_list1)
    child_mol2_dum = add_subs2(blank_mol, child_full_list2)
    child_mol1_nring_dum = dummy_atom_remover_non_ring(child_mol1_dum)
    child_mol2_nring_dum = dummy_atom_remover_non_ring(child_mol2_dum)
    child_mol1 = dummy_atom_remover_ring(child_mol1_nring_dum)
    child_mol2 = dummy_atom_remover_ring(child_mol2_nring_dum)
    children.append(child_mol1)
    children.append(child_mol2)
    return children

def void_uniform(mol1,mol2):
    import random
    children = []
    mol1_subs_list = get_subs(mol1)
    mol2_subs_list = get_subs(mol2)
    child_full_list1 = []
    child_full_list2 = []
    for x in range(0,10):
        p = random.randint(0,100)
        #print p
        if p <= 50:
            child_full_list1.append(mol1_subs_list[x])
            child_full_list2.append(mol2_subs_list[x])
        if p > 50:
            child_full_list1.append(mol2_subs_list[x])
            child_full_list2.append(mol1_subs_list[x])
    blank_mol = Chem.MolFromSmiles('Clc1cc(Cl)c(-c2cc(Cl)c(Cl)cc2Cl)cc1Br')
    child_mol1_dum = add_subs2(blank_mol, child_full_list1)
    child_mol2_dum = add_subs2(blank_mol, child_full_list2)
    child_mol1_nring_dum = dummy_atom_remover_non_ring(child_mol1_dum)
    child_mol2_nring_dum = dummy_atom_remover_non_ring(child_mol2_dum)
    child_mol1 = dummy_atom_remover_ring(child_mol1_nring_dum)
    child_mol2 = dummy_atom_remover_ring(child_mol2_nring_dum)
    children.append(child_mol1)
    children.append(child_mol2)
    return children

def void_crossover(p1,p2):
    import random
    p_smiles1 = Chem.MolToSmiles(p1)
    p_smiles2 = Chem.MolToSmiles(p2)
    print p_smiles1+" crossed with "+p_smiles2
    child_mols = []
    crossover_choice = random.randint(0,100)
    print crossover_choice
    if crossover_choice <= 30:
        print "crossover is single point"
        child_mols = void_single_point(p1,p2)
    if crossover_choice > 30 and crossover_choice < 90:
        print "crossover is double point"
        child_mols = void_double_point(p1,p2)
    if crossover_choice >= 90:
        print "crossover is uniform"
        child_mols = void_uniform(p1,p2)
    kids = []
    for c in child_mols:
        mutation = random.randint(0,100)
        if mutation >= 95:
            print "mutating"
            try:
                mutant_mol = void_point_mutation(c)
                kids.append(mutant_mol)
            except RuntimeError:
                print "mutation failed"
                kids.append(c)
        else:
            kids.append(c)
    return kids

def void_point_mutation(mol):
    import random
    subs_list = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,7,7,7,6,6,6,6,6,6,6,9,17,35,53]
    carbon_subs_smiles = ['C[Si](C)(C)(C)','C(C)(C)(C)','C[Si](CC)(CC)(CC)','C(F)(F)(F)'
                      ,'C#N','C#C', 'C']
    nitrogen_subs_smiles = ['*[H]','*(C)(C)','*(=O)(=O)']
    molh = Chem.AddHs(mol)
    rings = molh.GetRingInfo()
    ring_ats = rings.AtomRings()
    flat_ring_ats =  sum(ring_ats, ())
    list_atom_sites = []
    for ring_at in flat_ring_ats:
        n = [n.GetIdx() for n in molh.GetAtomWithIdx(ring_at).GetNeighbors()]
        non_ring_n = [nrn for nrn in n if nrn not in flat_ring_ats]
        list_atom_sites.append(non_ring_n)
    actual_atom_sites = [atom for sublist in list_atom_sites for atom in sublist]
    at_s = random.choice(actual_atom_sites)
    replacement = random.choice(subs_list)
    print at_s, replacement
    n = [n.GetIdx() for n in molh.GetAtomWithIdx(at_s).GetNeighbors()]
    print n
    molhe = Chem.EditableMol(mol)
    molhe.ReplaceAtom(at_s,Chem.Atom(0))
    if len(n) > 1:
        for i,neigh in enumerate(n):
            if i == 0:
                pass
            else:
                molhe.RemoveBond(at_s,neigh)
    temp = molhe.GetMol()
    Chem.SanitizeMol(temp)
    temp2 = Chem.GetMolFrags(temp, asMols=True)[0]
    molh = temp2
    if replacement == 6:
        print "carbon"
        carbon_sub = random.choice(carbon_subs_smiles)
        if carbon_sub == 'C':
            pmHe = Chem.EditableMol(molh)
            pmHe.ReplaceAtom(at_s, Chem.Atom(6))
        else:
            temp_mol = Chem.MolFromSmiles(carbon_sub)
            comb_mol = Chem.CombineMols(molh,temp_mol)
            frags = Chem.GetMolFrags(comb_mol)
            pmHe = Chem.EditableMol(comb_mol)
            pmHe.ReplaceAtom(at_s, Chem.Atom(6))
            if '[Si]' in carbon_sub:
                pmHe.AddBond(at_s,frags[1][0],Chem.BondType.TRIPLE)
            elif '#' in carbon_sub:
                pmHe.AddBond(at_s,frags[1][1],Chem.BondType.TRIPLE)
                pmHe.RemoveAtom(frags[1][0])
            else:
                pmHe.AddBond(at_s,frags[1][1],Chem.BondType.SINGLE)
                pmHe.AddBond(at_s,frags[1][2],Chem.BondType.SINGLE)
                pmHe.AddBond(at_s,frags[1][3],Chem.BondType.SINGLE)
                pmHe.RemoveAtom(frags[1][0])
    if replacement == 7:
        print "nitrogen"
        nitro_sub = random.choice(nitrogen_subs_smiles)
        print nitro_sub
        temp_mol = Chem.MolFromSmiles(nitro_sub)
        comb_mol = Chem.CombineMols(molh,temp_mol)
        frags = Chem.GetMolFrags(comb_mol)
        pmHe = Chem.EditableMol(comb_mol)
        pmHe.ReplaceAtom(at_s,Chem.Atom(7))
        if "O" in nitro_sub:
            pmHe.AddBond(at_s,frags[1][1],Chem.BondType.DOUBLE)
            pmHe.AddBond(at_s,frags[1][2],Chem.BondType.DOUBLE)
            pmHe.RemoveAtom(frags[1][0])
        elif "C" in nitro_sub:
            pmHe.RemoveAtom(frags[1][0])
            pmHe.AddBond(at_s,frags[1][0],Chem.BondType.SINGLE)
            pmHe.AddBond(at_s,frags[1][1],Chem.BondType.SINGLE)
        else:
            pmHe.ReplaceAtom(at_s,Chem.Atom(7))
            pmHe.RemoveAtom(frags[1][0])
    else:
        pmHe = Chem.EditableMol(molh)
        pmHe.ReplaceAtom(at_s,Chem.Atom(replacement))   
    tm2 = pmHe.GetMol()
    Chem.SanitizeMol(tm2)
    return tm2    

