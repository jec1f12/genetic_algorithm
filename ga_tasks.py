##this is the task file that is also replicated on iridis, it's here in case a run is ever performed on the same computer most of the functions are self explanatory
def com_writer(com_file,name):
    '''com file is an object containing all the coords and gauss arguments'''
    com = name+".com"
    com_f = open(com, "w")
    for line in com_file:
        if str(line).startswith("("):
             com_f.write("%6s %13.6f %13.6f %13.6f\n"% (line[0],float(line[1][0]),float(line[1][1]),float(line[1][2])))
        else:
             com_f.write(str(line))
    com_f.write(3*"\n")
    com_f.close()
    return com 

def energy_grabber(outfile):
    '''grabs the final energy from a gauss log file'''
    f = open(outfile, "r")
    f_list = list(f)
    f_list.reverse()
    i = 0
    for line in f_list:
        if "SCF Done" in line:
            energy = line.strip().split()[4]
            i += 1
        if i > 0:
            break
    return energy

def structure_grabber(outfile):
    '''grabs an optimised structure from a gauss log file'''
    f = open(outfile, "r")
    f_list = list(f)
    f_list.reverse()
    save_line = False
    struct_list = []
    for line in f_list:
        if "Input orientation" in line:
            break
        if "Distance matrix (angstroms)" in line:
            save_line = True
        if save_line == True:
            l = line.strip().split()
            struct_list.append(l)
    coords_list = []
    atom_no_list = []
    for l in struct_list:
        if len(l) < 5:
            pass
        try:
            int(l[0])
        except ValueError:
            pass
        else:
            coords_list.append(l[3:])
            atom_no_list.append(l[1])
    atom_list = []
    for atn in atom_no_list:
        if atn == "1":
            atom_list.append("H")
        if atn == "7":
            atom_list.append("N")
        if atn == "6":
            atom_list.append("C")
    structure = [(a1, c1) for a1,c1 in zip(atom_list,coords_list)]
    structure.reverse()
    return structure
 
def hl_grabber(outfile):
    '''gets homo lumo difference'''
    f = open(outfile, "r")
    f_list = list(f)
    f_list.reverse()
    save_line = False
    i = 0
    for j,line in enumerate(f_list):
        if "Population analysis using the SCF density." in line:
            break
        if "Alpha  occ. eigenvalues" in line:
            homo = line.split()[-1]
            lumo = f_list[j-1].split()[4]
            i += 1
        if i > 0:
            break
    hl_diff = float(homo) - float(lumo)
    return hl_diff

def lumo_grabber(outfile):
    '''gets lumo energy'''
    f = open(outfile, "r")
    f_list = list(f)
    f_list.reverse()
    save_line = False
    i = 0
    for j,line in enumerate(f_list):
        if "Population analysis using the SCF density." in line:
            break
        if "Alpha  occ. eigenvalues" in line:
            lumo = f_list[j-1].split()[4]
            i += 1
        if i > 0:
            break
    lumo_en = lumo
    return lumo_en

def max_d_grabber(outfile):
    '''gets maximum dipole moment'''
    f = open(outfile, "r")
    f_list = list(f)
    f_list.reverse()
    save_line = False
    i = 0
    for j,line in enumerate(f_list):
        if "Dipole moment (field-independent basis, Debye)" in line:
            max_dip = f_list[j-1].strip().split()[-1]
            i += 1
        if i > 0:
            break
    max_d = max_dip
    return max_d

