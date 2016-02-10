from celery import Celery
import subprocess as sb
import os
import shutil
import ga_tasks


app = Celery('tasks', backend='amqp' ,broker='amqp://myuser:mypassword@localhost:20100/myvhost')
import celeryconfig
app.config_from_object(celeryconfig)


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

