import os
import shutil
import signal
import celery
from celery import group
from tasks import *
import ga_setup_selection
import csp
import subprocess as sb

def mult_generation(all_mol_coords):
    '''generates multipole files for all population members'''
    xyz_files = ga_setup_selection.xyz_generator(all_mol_coords)
    job = group([run_mult.s(xyz) for xyz in xyz_files])#run_gauss is a celery task
    run_job = job.apply_async()#this is what runs the job of all the com files, none are returned till they're all done
    result1 = run_job.get()#this returns the results
    return result1

def mol_csp(result):
    complete = []
    for r in result:
        name = r[0]
        print "name is "+name
        isfilethere = r[1]
        print "is file there: "+str(isfilethere)
        xyz = r[2]
        print "xyz is "+str(xyz)
        if isfilethere == True:
            try:
                os.mkdir(name+"_server")
            except OSError:
                shutil.rmtree(name+"_server")
                os.mkdir(name+"_server")
            os.chdir(name+"_server")
            xyz_f = ga_setup_selection.res_writer(xyz,name)
            csp_command = """
                                python /home/jec1f12/CSPy_test/cspy/cspy/csp.py -s -i """+xyz_f+""" -exp -sat --portnum 20101 --server_version 1 -mvp 2.5 -sg 14 -zp 1 \
                                -ns 5 -mf ../"""+name+"""_mp.smult --dont_check_anisotropy t -pf /scratch/jec1f12/CSPy_bt/CSPy/potentials/w99rev_6311_s.pots \
                                --remove_negative_eigens --wcal_ngcv y --vdw_cutoff 25.0 --dmacrys_timeout 600.0 -zgs > server.out 2> server.err
                               """
            p = sb.Popen((csp_command),shell=True,  preexec_fn=os.setsid)
            print "server started"
            #p.communicate()
            job = run_csp.s(name)
            run_job = job.apply_async()
            result1 = run_job.get()
            complete.append(result1)
            os.killpg(os.getpgid(p.pid), signal.SIGTERM)
        os.chdir("..")
    return complete


def add_mol_csp(result,gm_names):
    complete = []
    for i,r in enumerate(result):
        name = r[0]
        print "name is "+name
        gm_name = gm_names[i]
        print "highest sobol name is "+gm_name
        isfilethere = r[1]
        print "is file there: "+str(isfilethere)
        xyz = r[2]
        print "xyz is "+str(xyz)
        prev_sobol_seed = gm_name.split("_")[4]
        print "previous highest sobol seed "+prev_sobol_seed
        sobol_seed = str(int(prev_sobol_seed) + 1)
        print "sobol seed for restart is "+sobol_seed
        if isfilethere == True:
            try:
                os.mkdir(name+"_server")
            except OSError:
                shutil.rmtree(name+"_server")
                os.mkdir(name+"_server")
            os.chdir(name+"_server")
            xyz_f = ga_setup_selection.res_writer(xyz,name)
            csp_command = """
                                python /home/jec1f12/CSPy_test/cspy/cspy/csp.py -s -i """+xyz_f+""" -exp -sat --portnum 20101 --server_version 1 -mvp 2.5 -sg 14 -zp 1 \
                                -ns 5 -mf ../"""+name+"""_mp.smult -ss """+sobol_seed+""" --dont_check_anisotropy t -pf /scratch/jec1f12/CSPy_bt/CSPy/potentials/w99rev_6311_s.pots \
                                --remove_negative_eigens --wcal_ngcv y --vdw_cutoff 25.0 --dmacrys_timeout 600.0 -zgs > server.out 2> server.err
                               """
            p = sb.Popen((csp_command),shell=True,  preexec_fn=os.setsid)
            print "server started"
            #p.communicate()
            job = run_csp.s(name)
            run_job = job.apply_async()
            result1 = run_job.get()
            complete.append(result1)
            os.killpg(os.getpgid(p.pid), signal.SIGTERM)
        os.chdir("..")
    return complete







 
