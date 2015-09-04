from celery import Celery
import celeryconfig
import subprocess as sb
import os
import shutil

app = Celery('tasks', backend='amqp' ,broker='amqp://myuser:mypassword@localhost:20100/myvhost')
app.config_from_object(celeryconfig)


@app.task
def add(x, y):
    return x + y


@app.task
def mul(x, y):
    return x * y


@app.task
def xsum(numbers):
    return sum(numbers)

@app.task
def run_gauss(com_file):
    name = com_file[5].split(" ")[-1].strip()
    com = name+".com"
    com_f = open(com, "w")
    for line in com_file:
        if str(line).startswith("("):
             com_f.write("%6s %13.6f %13.6f %13.6f\n"% (line[0],float(line[1][0]),float(line[1][1]),float(line[1][2])))
        else:
             com_f.write(str(line))
    com_f.write(3*"\n")
    com_f.close()
    run_dir = name+"_run"
    os.mkdir(run_dir)
    sb.call(("cp "+com+" "+run_dir),shell=True)
    os.chdir(run_dir)
    sb.call(("g09 "+com),shell=True)
    outfile = name+".log"
    f = open(outfile, "r")
    f_list = list(f)
    f_list.reverse()
    save_line = False
    i = 0
    for line in f_list:
        if "SCF Done" in line:
            energy = line.strip().split()[4]
            i =+ 1
        if i > 0:
            break
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
    #print coords_list
    atom_list = []
    for atn in atom_no_list:
        if atn == "1":
            atom_list.append("H")
        if atn == "7":
            atom_list.append("N")
        if atn == "6":
            atom_list.append("C")
    #print atom_list
    structure = [(a1, c1) for a1,c1 in zip(atom_list,coords_list)]
    structure.reverse()
    os.chdir("..")
    #shutil.rmtree(run_dir)
    return (name,energy,structure)

@app.task
def run_gauss2(com_file):
    name = com_file[5].split(" ")[-1].strip()
    new_name = name+"_opt"
    com = name+".com"
    com_f = open(com, "w")
    for line in com_file:
        if str(line).startswith("("):
             com_f.write("%6s %13.6f %13.6f %13.6f\n"% (line[0],float(line[1][0]),float(line[1][1]),float(line[1][2])))
        else:
             com_f.write(str(line))
    com_f.write(3*"\n")
    com_f.close()
    run_dir = name+"_run"
    os.mkdir(run_dir)
    sb.call(("cp "+com+" "+run_dir),shell=True)
    os.chdir(run_dir)
    sb.call(("g09 "+com),shell=True)
    outfile = name+".log"
    f = open(outfile, "r")
    f_list = list(f)
    f_list.reverse()
    save_line = False
    i = 0
    for line in f_list:
        if "SCF Done" in line:
            energy = line.strip().split()[4]
            i += 1
        if i > 0:
            break
    os.chdir("..")
    #shutil.rmtree(run_dir)
    return (name,energy) 
