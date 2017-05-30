import random
import copy
from ga_setup_selection import mutator
from rdkit import Chem

####atomic site crossover functions, this crossover uses positions on outside of fused ring system as the genome to be crossed over



def crossover_chooser(crossover_choice):
    '''just returns a value based on the random numbers for crossover_choice'''
    if crossover_choice <= 30:
        return 1
    if crossover_choice > 30:
        if crossover_choice < 90:
            return 2
    if crossover_choice >= 90:
        return 3



def smiles_fixer(smiles, bit_string, blank_pos):
    '''some smiles are written in an incompatible format for the crossover pair function
    usually happens when 2Ns are in postions 0,2 and 16,17'''
    bad_crossover = []
    bad_pos = [0,2,3,4,7,10,13,16,17,18,19,22,25,28]
    bad_smiles = smiles
    good_pos = bit_string
    blanks = blank_pos
    for i,atom in enumerate(bad_smiles):
        for bp in bad_pos:
            if i == bp:
                bad_crossover.append(atom)
    good_cross = [(i1, i2) for i1, i2 in zip(good_pos, bad_crossover)]
    new_mask = good_cross + blank_pos
    new_mask_sorted = sorted(new_mask, key = lambda x: int(x[0]))
    new_mask_list = [atom for pos,atom in new_mask_sorted]
    test_smile = ''.join(new_mask_list)
    return test_smile


def crossover_pair_varying(mol_1,mol_2):
    '''this function crosses over 2 seperate molecules, currently 3 types of crossover are
    implemented, single point, two point and uniform, which operator is applied is decided 
    randomly, this function will then return 2 child molecules'''
    children = [] #child molecules will be appended to this list
    bit_string = dict.fromkeys(['0','2','3','6','9','12','15','16','17','18','21','24',
                            '27','30']) #acceptable atom positions to crossover
    blank_pos = [('4','c'),('7','c'), ('10','c'),('13','c'),('19','c'),
             ('22','c'),('25','c'),('28','c'),('1','1'),('5','2'),('8','3'),
             ('11','4'),('14','5'),('20','5'),('23','4'),('26','3'),('29','2'),
             ('31','1')] #ring central atoms and numbers in the smiles string 
    standard_smiles = 'c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1' #pentacene smiles
    smiles_1 = Chem.MolToSmiles(mol_1) 
    smiles_2 = Chem.MolToSmiles(mol_2)
    #print other_bit_string
    if smiles_1.endswith('2'): #some smiles can break this function usually < 1%
        smiles_1 = smiles_fixer(smiles_1,bit_string,blank_pos)
    if smiles_2.endswith('2'):
        smiles_2 = smiles_fixer(smiles_2,bit_string,blank_pos)
    bit_string_1 = copy.deepcopy(bit_string) #stops me breaking the original
    bit_string_2 = copy.deepcopy(bit_string)
    smiles_list_1 = list(smiles_1) #usually easier to work with the list form of the smiles
    smiles_list_2 = list(smiles_2)
    stan_smiles_list = list(standard_smiles)
    #print bit_string_1
    for pos in bit_string_1: #this adds the correct atoms to the bit string for both smiles
        for i,atom in enumerate(smiles_list_1):
            if int(pos) == i:
                bit_string_1.update({pos:atom})
    for pos in bit_string_2:
        for i,atom in enumerate(smiles_list_2):
            if int(pos) == i:
                bit_string_2.update({pos:atom})
    bts_1_sorted = [] #now sort the dict by the position 
    for key in sorted(bit_string_1, key = int):
        bts_1_sorted.append((key, bit_string_1[key]))
    #print bts_1_sorted
    bts_2_sorted = []
    for key in sorted(bit_string_2, key = int):
        bts_2_sorted.append((key, bit_string_2[key]))
    #print bts_2_sorted
    mask = [] #this is what will contain our new molecules
    mask2 = []
    off_count = [1,2] #just values to make sure we make 2 unqiue offspring 
    crossover_choice = random.randint(0,100) #choosing the number that will decide the crossover type
    x = crossover_chooser(crossover_choice) #send it to the function above
    if x == 1:
        print "crossover is single_point"
        n = random.randint(1,7) #where we'll crossover the molecules
        print "point is "+str(n)
        for off in off_count:
            if off == 1:
                for i,atom in enumerate(stan_smiles_list): #reading along the full standard smiles 
                    for pos,atom2 in bts_1_sorted[:n]: #and the sorted bit string (0, c) etc upto our crossover point
                        if i == int(pos): #if the position in the standard smiles == the position of our crossover atoms
                            mask.append((pos, atom2)) #append the position and atom to our new molecule
                    for pos2,atom3 in bts_2_sorted[n:]: #reading from the second molecule from the crossover point
                        if i == int(pos):
                            mask.append((pos2,atom3)) #append position and atom
            #children.append(mask)
            if off == 2: #this is the reverse of the code above, ie start with parent 2
                for i,atom in enumerate(stan_smiles_list):
                    for pos3,atom4 in bts_2_sorted[:n]:
                        if i == int(pos):
                            mask2.append((pos3,atom4))
                    for pos4,atom5 in bts_1_sorted[n:]:
                        if i == int(pos):
                            mask2.append((pos4,atom5))
        children.append(mask) #add both masks to the child list 
        children.append(mask2) #at this point masks should be a full list of 14 atoms and positions     
    if  x == 2:
        print "crossover is two_point"
        n0 = random.randint(1,7) #two point crossover needs two points, n0 for parent 1, n1 for parent2
        n1 = random.randint(1,7)
        print "Taking "+str(n0)+" of parent 1"
        print "Taking "+str(n1)+" of parent 2"
        for off in off_count:
            if off == 1:
                for i,atom in enumerate(stan_smiles_list):
                    for pos,atom2 in bts_1_sorted[:n0]: #very similar to above take upto n0
                        if i == int(pos):
                            mask.append((pos,atom2))
                    for pos2,atom3 in bts_2_sorted[n0:(n0+n1)]: #differs here, we only need n1 length bit of parent 2
                        if i == int(pos):
                            mask.append((pos2,atom3))
                n3 = len(mask) #depending on values of n we can be left short a full mask
                if n3 < 14: 
                    for i,atom in enumerate(stan_smiles_list):
                        for pos,atom2 in bts_1_sorted[n3:]: #two point crossover needs the rest of parent 1, taken from end of mask
                            if i == int(pos):
                                mask.append((pos,atom2))
                children.append(mask)
            if off == 2: #the same as above just reversed parent order
                for i,atom in enumerate(stan_smiles_list):
                    for pos,atom2 in bts_2_sorted[:n0]:
                        if i == int(pos):
                            mask2.append((pos,atom2))
                    for pos2,atom3 in bts_1_sorted[n0:(n0+n1)]:
                        if i == int(pos):
                            mask2.append((pos2,atom3))
                n3 = len(mask2)
                if n3 < 14:
                    for i,atom in enumerate(stan_smiles_list):
                        for pos,atom2 in bts_2_sorted[n3:]:
                            if i == int(pos):
                                mask2.append((pos,atom2))
                children.append(mask2)
    if x == 3:
        print "crossover is uniform"
        pos_list = range(0,14) #uniform crossover is a random chance for each point being either from parent 1 or parent 2
        #print len(bts_1_sorted)
        for off in off_count:
            if off == 1:
                for i in pos_list:
                    #print i
                    p = random.randint(0,100)
                    if p > 50: #50% chance is the standard 
                        mask.append(bts_1_sorted[i]) #this ensures the 2 children will be mirror images of each other
                        mask2.append(bts_2_sorted[i])
                    if p <= 50:
                        mask2.append(bts_1_sorted[i])
                        mask.append(bts_2_sorted[i])
        children.append(mask)
        children.append(mask2)
            #print children
    child_mols = []
    for child in children:
        full_child_unsorted = blank_pos + child #create full list for smiles
        full_child = sorted(full_child_unsorted, key = lambda x: int(x[0])) #sort by position
        child_smiles_list = [atom for pos,atom in full_child] #get rid of position
        child_smile = ''.join(child_smiles_list) #turn into smiles
        child_mol = Chem.MolFromSmiles(child_smile) #turn into full molecule
        child_mols.append(child_mol)
    for cm in child_mols: #just error checking for now, any junk molecule will fail the turning into mol
        if cm == None:
            print smiles_1
            print mask
            #print mask2
            print bts_1_sorted
            #print smiles_2
            #print bts_2_sorted
    #print len(child_mols)
    return child_mols

def smiles_fixer_fixed(smiles):
    '''some smiles also break the fixed crossover, crossovers are too different to share the same fixer'''
    smiles_list = list(smiles)
    atoms = []
    number_pos = [('1', 0), ('2', 3), ('3',5), ('4',7),('5',9),
              ('5',14), ('4',16),('3',18),('2',20),('1',22)] 
    j = 0
    for i,atom in enumerate(smiles_list):
        if i >= 15:
            try:
                int(atom)
            except ValueError:
                atoms.append((atom,j))
                j += 1
    for i,atom in enumerate(smiles_list):
        if i < 15:
            try:
                int(atom)
            except ValueError:
                atoms.append((atom, j))
                j +=1
    new_unsorted = atoms+number_pos
    new_sorted = sorted(new_unsorted, key = lambda x: int(x[1]))
    new_smiles_list = [atom for atom,pos in new_sorted]
    new_smile = ''.join(new_smiles_list)
    return new_smile

def n_counter(mol):
    '''returns a list of N positions, len of list will be used for number checking'''
    n_pos = []
    smiles = Chem.MolToSmiles(mol)
    smiles_list = list(smiles)
    for i,atom in enumerate(smiles):
        if atom == 'n':
            n_pos.append(i)
    return n_pos

def crossover_pair_fixed(mol_1,mol_2):
    '''this function crosses over two molecules using only the N positions as the genome, allowing the N to be fixed
       usually GA is run with varying'''   
    n_number = 5
    parent_1 =  Chem.MolToSmiles(mol_1)
    parent_2 =  Chem.MolToSmiles(mol_2)
    print "parent 1 = "+parent_1
    print "parent_2 = "+parent_2
    blank = 'c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1'
    blank_list = list(blank)
    off_count = [1,2]
    if parent_1[-3] != '2':
        parent_1 = smiles_fixer_fixed(parent_1)
        print "new parent_1 "+parent_1
    if parent_2[-3] != '2':
        parent_2 = smiles_fixer_fixed(parent_2)
        print "new parent_2 "+parent_2
    p1_list = list(parent_1)
    p2_list = list(parent_2)
    p1_n = []
    p2_n = []
    for i,atom in enumerate(p1_list):
        if atom == 'n':
            p1_n.append(i)
    for i,atom in enumerate(p2_list):
        if atom == 'n':
            p2_n.append(i)
    print p1_n
    print p2_n
    crossover_choice = random.randint(0,100) #choosing the number that will decide the crossover type
    crossover_type =  crossover_chooser(crossover_choice) #send it to the function above
    off_count = [1,2]
    children = []
    mask_1 = []
    child_1 = []
    mask_2 = []
    child_2 = []
    if crossover_type == 1:    
        print "crossover is single point"
        n = random.randint(1,n_number)
        print n
        for off in off_count:
            if off == 1:
                for pos in p1_n[:n]:
                    mask_1.append(pos)
                for pos in p2_n[n:]:
                    mask_1.append(pos)
                #print mask_1
                for i,atom in enumerate(blank_list):
                    if i in mask_1:
                        atom = 'n'
                    child_1.append(atom)
                #print child_1
            if off == 2:
                for pos in p2_n[:n]:
                    mask_2.append(pos)
                for pos in p1_n[n:]:
                    mask_2.append(pos)
                #print mask_2
                for i,atom in enumerate(blank_list):
                    if i in mask_2:
                        atom = 'n'
                    child_2.append(atom)
                #print child_2
    if crossover_type == 2:
        print "crossover is double point"
        n0 = random.randint(1,n_number)
        n1 = random.randint(1,n_number)
        print n0
        print n1
        for off in off_count:
            if off == 1:
                for pos in p1_n[:n0]:
                    mask_1.append(pos)
                for pos in p2_n[n0:(n0+n1)]:
                    mask_1.append(pos)
                n3 = len(mask_1)
                if n3 < n_number:
                    for pos in p1_n[n3:]:
                        mask_1.append(pos)
                #print mask_1
                for i,atom in enumerate(blank_list):
                    if i in mask_1:
                        atom = 'n'
                    child_1.append(atom)
                #print child_1
            if off == 2:
                for pos in p2_n[:n0]:
                    mask_2.append(pos)
                for pos in p1_n[n0:(n0+n1)]:
                    mask_2.append(pos)
                n3 = len(mask_2)
                if n3 < n_number:
                    for pos in p2_n[n3:]:
                        mask_2.append(pos)
                #print mask_2
                for i,atom in enumerate(blank_list):
                    if i in mask_2:
                        atom = 'n'
                    child_2.append(atom)
                #print child_2
    if crossover_type == 3:
        print "crossover is uniform"
        pos_list = range(0,5)
        for off in off_count:
            if off == 1:
                for i in pos_list:
                    p = random.randint(0,100)
                    if p >= 50:
                        mask_1.append(p1_n[i])
                        mask_2.append(p2_n[i])
                    if p < 50:
                        mask_2.append(p1_n[i])
                        mask_1.append(p2_n[i])
                #print mask_1
                #print mask_2
        for i,atom in enumerate(blank_list):
            if i in mask_1:
                atom = 'n'
            child_1.append(atom)
        for i,atom in enumerate(blank_list):
            if i in mask_2:
                atom = 'n'
            child_2.append(atom)            
    children.append(child_1)
    children.append(child_2)
    child_mols = []
    for child in children:
        child_smiles = ''.join(child)
        print child_smiles
        ch_mol = Chem.MolFromSmiles(child_smiles)
        child_mols.append(ch_mol)
        if ch_mol == None:
            print parent_1+" + "+parent_2
    
    return child_mols
