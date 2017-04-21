from celery import Celery
import subprocess as sb
import os
import shutil
import ga_tasks


app = Celery('tasks', backend='amqp' ,broker='amqp://myuser:mypassword@localhost:20100/myvhost3')
import celeryconfig
app.config_from_object(celeryconfig)


@app.task
def run_mult(xyz_file):
    print xyz_file
    name = xyz_file[0]
    xyz_coords = xyz_file[1]
    run_dir = name+"_run"
    try:
        os.mkdir(run_dir)
    except OSError:
        shutil.rmtree(run_dir)
        os.mkdir(run_dir)
    xyz = ga_tasks.xyz_writer(xyz_coords,name)
    sb.call(("mv "+xyz+" "+run_dir),shell=True)
    os.chdir(run_dir)
    print os.getcwd()
    sb.call(("bash /scratch/jec1f12/CSPy_bt/CSPy/pbs/mult.sh"),shell=True)
    isfilethere = os.path.isfile(name+"_mp.smult")
    res_file = ga_tasks.res_grabber(name)
    os.chdir("..")
    #shutil.rmtree(run_dir)
    return (name,isfilethere,res_file)


@app.task
def run_csp(result):
    print name
    print name+" csp started"
    print os.getcwd()
    run_dir = name+"_test_run"
    os.chdir(run_dir)
    print os.getcwd()
    #connect_cmd = "ssh -N -L 19753:localhost:20101 152.78.197.25 &"
    connect = sb.Popen(("ssh -N -L 19753:localhost:20101 152.78.197.25"),shell=True, preexec_fn=os.setsid)
    print "connection made"
    py_cmd = "python /scratch/jec1f12/CSPy_bt/CSPy/cspy/csp.py -c -ncp 12 -ip 127.0.0.1 -wt 02:00:00  > worker_1.out 2> worker_1.err"
    sb.call((py_cmd),shell=True)
    print "command run"
    ga_tasks.move_results(name)
    global_minimum = ga_tasks.get_global_minimum(name)
    gm_name = ga_tasks.get_gm_name(name)
    os.killpg(os.getpgid(connect.pid), signal.SIGTERM)
    os.chdir("..")
    shutil.rmtree(name+"_run")
    print global_minimum
    return (name,global_minimum,gm_name)




@app.task
def run_gauss(com_file):
    name = com_file[5].split(" ")[-1].strip()
    run_dir = name+"_run"
    try:
        os.mkdir(run_dir) 
    except OSError:
        shutil.rmtree(run_dir)
        os.mkdir(run_dir)    
    com = ga_tasks.com_writer(com_file,name)
    sb.call(("mv "+com+" "+run_dir),shell=True)
    os.chdir(run_dir)
    sb.call(("g09 "+com),shell=True)
    outfile = name+".log"
    energy = ga_tasks.energy_grabber(outfile)
    structure = ga_tasks.structure_grabber(outfile)
    os.chdir("..")
    shutil.rmtree(run_dir)
    return (name,energy,structure)

@app.task
def run_gauss2(com_file):
    name = com_file[5].split(" ")[-1].strip()
    run_dir = name+"_run"
    try:
        os.mkdir(run_dir)
    except OSError:
        shutil.rmtree(run_dir)
        os.mkdir(run_dir)
    com = ga_tasks.com_writer(com_file,name)
    sb.call(("mv "+com+" "+run_dir),shell=True)
    os.chdir(run_dir)
    sb.call(("g09 "+com),shell=True)
    outfile = name+".log"
    energy = ga_tasks.energy_grabber(outfile)
    os.chdir("..")
    shutil.rmtree(run_dir)
    return (name,energy) 

@app.task
def run_gausshl(com_file):
    name = com_file[5].split(" ")[-1].strip()
    run_dir = name+"_run"
    try:
        os.mkdir(run_dir)
    except OSError:
        shutil.rmtree(run_dir)
        os.mkdir(run_dir)
    com = ga_tasks.com_writer(com_file,name)
    sb.call(("mv "+com+" "+run_dir),shell=True)
    os.chdir(run_dir)
    sb.call(("g09run "+com),shell=True)
    outfile = name+".log"
    hl_diff = ga_tasks.hl_grabber(outfile)
    os.chdir("..")
    shutil.rmtree(run_dir)
    return(name,hl_diff)

@app.task
def run_gaussl(com_file):
    name = com_file[5].split(" ")[-1].strip()
    run_dir = name+"_run"
    try:
        os.mkdir(run_dir)
    except OSError:
        shutil.rmtree(run_dir)
        os.mkdir(run_dir)
    com = ga_tasks.com_writer(com_file,name)
    sb.call(("mv "+com+" "+run_dir),shell=True)
    os.chdir(run_dir)
    sb.call(("g09run "+com),shell=True)
    outfile = name+".log"
    hl_diff = ga_tasks.lumo_grabber(outfile)
    os.chdir("..")
    shutil.rmtree(run_dir)
    return(name,lumo_en)

@app.task
def run_gaussmd(com_file):
    name = com_file[5].split(" ")[-1].strip()
    run_dir = name+"_run"
    try:
        os.mkdir(run_dir)
    except OSError:
        shutil.rmtree(run_dir)
        os.mkdir(run_dir)
    com = ga_tasks.com_writer(com_file,name)
    sb.call(("mv "+com+" "+run_dir),shell=True)
    os.chdir(run_dir)
    sb.call(("g09run "+com),shell=True)
    outfile = name+".log"
    max_d_moment = ga_tasks.max_d_grabber(outfile)
    os.chdir("..")
    shutil.rmtree(run_dir)
    return(name,max_d_moment)

