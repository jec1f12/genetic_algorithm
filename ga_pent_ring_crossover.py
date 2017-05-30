import random
from rdkit import Chem
from rdkit.Chem import AllChem
import random
def find_dummy_atoms(mol):
    '''finds the dummy atoms in a mol, returns atom ids'''
    dummy_atoms = []
    for atom in mol.GetAtoms():
        at_num = atom.GetAtomicNum()
        #print at_num
        if at_num == 0:
            index = atom.GetIdx()
            neigh = []
            neighbours = atom.GetNeighbors()
            for n in neighbours:
                neigh_id = n.GetIdx()
                neigh.append(neigh_id)
            dummy_atoms.append((index, neigh))
            
    return dummy_atoms

def fragment_total(mol):
    '''fragments molecule into all constituent rings, returns a list of each ring'''
    frags = []
    origmol = mol
    for x in range(0,5):#loop over number of rings
        if x == 4:
            tempmol = origmol
            cut_pts = return_cut_points(tempmol)#get all cutpoints in a molecule
            print cut_pts
            cp = get_actual_cut_pts(tempmol,(x),cut_pts)#get correct cutpoints for ring number
            print cp
            frag = frag_mol(tempmol,cp,bf='n')#break ring off molecule
            frags.append(frag[1])
            break
        if x == 0:#skip breaking first ring for now
            tempmol = origmol
        r = tempmol.GetRingInfo() #get first ring
        #print r.NumRings()
        cut_pts = return_cut_points(tempmol)
        print cut_pts
        cp = get_actual_cut_pts(tempmol,1,cut_pts)
        print cp
        frag = frag_mol(tempmol,cp,bf='n')
        frags.append(frag[0])
        tempmol = frag[1]
    print len(frags)
    return frags

def join_total(frags):
    '''joins list of multiple fragments together, returns new mol'''
    print len(frags)
    join_pts = [0,1,2,3] #similar loop to fragment total
    for x in join_pts:
        #print x
        if x == 0:
            temp_mol = Chem.CombineMols(frags[x],frags[x+1])#makes one molecule object out of two separate ones
            atch_mol = join_mols(temp_mol)
            #r = atch_mol.GetRingInfo()
            #print r.NumRings()
            continue
        if x == 4:
            break
        else:
            temp_mol = Chem.CombineMols(atch_mol, frags[x+1])
            new_atch_mol = join_mols(temp_mol)
            #r = new_atch_mol.GetRingInfo()
            #print r.NumRings()
            atch_mol = new_atch_mol
    return atch_mol



def ring_fixer(cross_mol,parent_mol):
    '''two point crossover can lead to too short and too long molecules,
    this function either adds or subtracts from the molecule'''
    rings = cross_mol.GetRingInfo()
    num_rings = rings.NumRings()
    print "number of rings = "+str(num_rings)
    if num_rings == 5:
        new_mol = cross_mol#molecule doesn't need fixing
    if num_rings < 5:
        ring_ats = rings.AtomRings()[num_rings-1]#get ring atoms for last ring
        end_ats = []
        for atom in cross_mol.GetAtoms():
            if atom.GetIdx() in ring_ats:
                if len(atom.GetBonds()) == 2:#get end atoms of last ring
                     end_ats.append(atom.GetIdx())
        cp = return_cut_points(parent_mol)
        ring_cp = [cp[num_rings-1] + cp[-num_rings]]
        print ring_cp
        actual_cp = get_actual_cut_pts(parent_mol,num_rings,cp)
        print actual_cp
        additional_ring = frag_mol(parent_mol,actual_cp,bf='n')[1]#get additional ring from parent mol
        te2_m = Chem.EditableMol(cross_mol)
        te2_m.ReplaceAtom(end_ats[1], Chem.Atom(0))#replacing with dummy atoms allows for joining
        te2_m.ReplaceAtom(end_ats[2], Chem.Atom(0))
        tm2 = te2_m.GetMol()
        Chem.SanitizeMol(tm2)
        temp_mol = Chem.CombineMols(tm2,additional_ring)
        new_mol = join_mols(temp_mol)
    if num_rings > 5:
        cut_pts = return_cut_points(cross_mol)
        ring_ats = rings.AtomRings()[5]#get ring ats for 5th ring
        cut_atoms = []
        cut_atoms.append(ring_ats[-1])
        cut_atoms.append(ring_ats[-2])
        print cut_atoms
        cp_tup = [(cut_at,cut_pt) for cut_at,cut_pt in cut_pts if cut_at in cut_atoms]
        frags = frag_mol(cross_mol,cp_tup,bf='n')#break at 5th ring
        temp_mol = frags[0]
        dum_ats = find_dummy_atoms(temp_mol)
        te_m = Chem.EditableMol(temp_mol)
        rem_at = te_m.RemoveAtom(dum_ats[1][0])#remove dummy atoms
        rem_at2 = te_m.RemoveAtom(dum_ats[0][0])
        new_at = te_m.AddAtom(Chem.Atom(6))## these new atoms are added onto the end of the atom idx
        new_at2 = te_m.AddAtom(Chem.Atom(6))##all orig atoms are moved down by 1
        te_m.AddBond(new_at,dum_ats[0][1][0],Chem.BondType.AROMATIC)
        te_m.AddBond(new_at,new_at2,Chem.BondType.AROMATIC)#need to remake aromatic bonds for new carbon atoms
        te_m.AddBond(new_at2,dum_ats[-1][-1][1],Chem.BondType.AROMATIC)
        tm = te_m.GetMol()
        Chem.SanitizeMol(tm)
        new_mol = tm
    return new_mol



def return_cut_points(mol):
    '''returns cutting points for crossover, gets all 3 center C 
    and returns specific ones for different rings'''
    cut_pts = []
    for atom in mol.GetAtoms(): ##can functionalise this snippet just to return cut_pts
        index = atom.GetIdx()
        at_num = atom.GetAtomicNum()
        neigh = atom.GetNeighbors()
        neighbours = []
        if len(neigh) > 2:
            for n in neigh:
                neigh_id = n.GetIdx()
                neighbours.append(neigh_id)
            cut_pts.append((index, neighbours))
    #actual_cp = cut_pts[num-1] + cut_pts[-num]
    return cut_pts

def frag_mol(mol,actual_cp,bf):
    '''frag mol based on cutting points, add dummy atoms at points'''
    frags = []
    print actual_cp
    atom1 = actual_cp[0][0]
    atom2 = actual_cp[1][0]
    atom_diff = abs(atom1 - atom2)#sometimes atom ordering is a bit iffy so atom diff controls bonds being created correctly
    #print atom_diff
    em = Chem.EditableMol(mol)
    if atom_diff == 1:
        em.RemoveBond(atom1,actual_cp[0][1][0])
        em.RemoveBond(atom2,actual_cp[1][1][1])
    else:
        em.RemoveBond(atom1,actual_cp[0][1][1])
        em.RemoveBond(atom2,actual_cp[1][1][0])
    new_at = em.AddAtom(Chem.Atom(0))
    if atom_diff == 1:
        em.AddBond(actual_cp[0][1][0],new_at,Chem.BondType.AROMATIC)
    else:
        em.AddBond(actual_cp[0][1][1],new_at,Chem.BondType.AROMATIC)
    new_at2 = em.AddAtom(Chem.Atom(0))
    if atom_diff == 1:
        em.AddBond(actual_cp[1][1][1],new_at2,Chem.BondType.AROMATIC)
    else:
        em.AddBond(actual_cp[1][1][0],new_at2,Chem.BondType.AROMATIC)
    em.AddBond(new_at,new_at2,Chem.BondType.AROMATIC)
    if bf == "y":#branch factor
        side = random.randint(1,2)
        print side#side controls if dummy atoms are moved up or down the ring
        if side == 1:
            em.ReplaceAtom(atom1, Chem.Atom(0))
            em.ReplaceAtom(actual_cp[0][1][0], Chem.Atom(0))
        if side == 2:
            em.ReplaceAtom(actual_cp[1][1][1], Chem.Atom(0))
            em.ReplaceAtom(atom2, Chem.Atom(0))
    if bf == 'n':
        em.ReplaceAtom(atom1, Chem.Atom(0))
        em.ReplaceAtom(atom2, Chem.Atom(0))
    nm = em.GetMol()
    frags = Chem.GetMolFrags(nm,asMols=True)
    print len(frags)
    return frags

def join_mols(mol):
    '''join fragged mols, based on dummy atom positions'''
    dummy_atoms = find_dummy_atoms(mol)
    #print dummy_atoms 
    #print dummy_atoms[0][0]
    temp_em = Chem.EditableMol(mol)
    temp_em.RemoveAtom(dummy_atoms[-1][0])
    temp_em.RemoveAtom(dummy_atoms[-2][0])
    temp_em.ReplaceAtom(dummy_atoms[1][0],Chem.Atom(6))
    temp_em.ReplaceAtom(dummy_atoms[0][0],Chem.Atom(6))
    temp_em.AddBond(dummy_atoms[0][0],dummy_atoms[-2][1][0],Chem.BondType.AROMATIC)
    temp_em.AddBond(dummy_atoms[1][0],dummy_atoms[-1][1][0],Chem.BondType.AROMATIC)
    tm = temp_em.GetMol()
    Chem.SanitizeMol(tm)
    return tm


def get_actual_cut_pts(mol,n,cut_pts):
    '''returns cut points for a specified ring'''
    rings = mol.GetRingInfo()
    ring_ats = rings.AtomRings()[n]
    actual_ca = ring_ats[-2:]
    actual_cp = [(at,neigh) for at,neigh in cut_pts if at in actual_ca]
    return actual_cp

def single_point_crossover(m1,m2,bf):
    '''single point crossover for rings'''
    crossed_mols = []
    n1 =  random.randint(1,4)
    cp1 = return_cut_points(m1)
    cp2 = return_cut_points(m2)
    actual_cp1 = get_actual_cut_pts(m1,n1,cp1) 
    actual_cp2 = get_actual_cut_pts(m2,n1,cp2)
    frags1 = frag_mol(m1,actual_cp1,bf)
    frags2 = frag_mol(m2,actual_cp2,bf)
    temp_mol1 = Chem.CombineMols(frags1[0],frags2[1])
    temp_mol2 = Chem.CombineMols(frags2[0],frags1[1])
    crossed_mol1 = join_mols(temp_mol1)
    crossed_mol2 = join_mols(temp_mol2)
    crossed_mols.append(crossed_mol1)
    crossed_mols.append(crossed_mol2)
    return crossed_mols

def double_point_crossover(m1,m2,bf):
    '''double point crossover for rings'''
    crossed_mols = []
    n1 =  random.randint(1,4)
    n2 =  random.randint(n1-1,4)
    print n1
    print n2
    if n2 == 0:
        n2 = 1
    cp1 = return_cut_points(m1)
    cp2 = return_cut_points(m2)
    actual_cp1 = get_actual_cut_pts(m1,n1,cp1)
    #print actual_cp1
    actual_cp2 = get_actual_cut_pts(m2,n2,cp2)
    #print actual_cp2
    frags1 = frag_mol(m1,actual_cp1,bf)
    frags2 = frag_mol(m2,actual_cp2,bf)
    print frags1
    print frags2
    fm =  [Chem.MolToSmiles(x,True) for x in frags1]
    fm2 = [Chem.MolToSmiles(x,True) for x in frags2]
    print fm
    print fm2
    temp_mol1 = Chem.CombineMols(frags1[0],frags2[1])
    temp_mol2 = Chem.CombineMols(frags2[0],frags1[1])
    crossed_mol1 = join_mols(temp_mol1)
    crossed_mol2 = join_mols(temp_mol2)
    fixed_mol1 = ring_fixer(crossed_mol1,m1)
    fixed_mol2 = ring_fixer(crossed_mol2,m2)
    crossed_mols.append(fixed_mol1)
    crossed_mols.append(fixed_mol2)
    return crossed_mols

def uniform_crossover(m1,m2):
    '''uniform crossover for rings'''
    crossed_mols = []
    new_mol1 = []
    new_mol2 = []
    frags1 = fragment_total(m1)
    frags2 = fragment_total(m2)
    #print frags1
    #print frags2
    for x in range(0,5):
        p = random.randint(0,100)
        #print p
        if p <= 50:
            new_mol1.append(frags1[x])
            new_mol2.append(frags2[x])
        if p > 50:
            new_mol1.append(frags2[x])
            new_mol2.append(frags1[x])
    crossed_mol1 = join_total(new_mol1)
    crossed_mol2 = join_total(new_mol2)
    crossed_mols.append(crossed_mol1)
    crossed_mols.append(crossed_mol2)
    return crossed_mols

def ring_crossover(p1,p2):
    '''crosses over 2 seperate molecules, using rings as the genome, also can create branched mols'''
    p_smiles1 = Chem.MolToSmiles(p1)
    p_smiles2 = Chem.MolToSmiles(p2)
    m1 = Chem.MolFromSmiles(p_smiles1)
    m2 = Chem.MolFromSmiles(p_smiles2)
    print p_smiles1+" crossed with "+p_smiles2
    child_mols = []
    branch_factor = random.randint(0,100)
    if branch_factor <= 5:
        bf = 'y'
    else:
        bf = 'n'
    if bf == 'y':
        print "branching"
    crossover_choice = random.randint(0,100)
    if "(" in p_smiles1 or p_smiles2:
        crossover_choice = random.randint(0,89)
    print crossover_choice 
    if crossover_choice <= 30:
        print "crossover is single point"
        child_mols = single_point_crossover(m1,m2,bf)
    if crossover_choice > 30 and crossover_choice < 90:
        print "crossover is double point"
        child_mols = double_point_crossover(m1,m2,bf)
    if crossover_choice >= 90:
        print "crossover is uniform"
        child_mols = uniform_crossover(m1,m2)
    return child_mols
