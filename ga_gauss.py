import celery
from celery import group
from tasks import *
import ga_setup_selection




def reorg_en_calculator(energies):
    '''calculate reorganisation energies from run results'''
    result = []
    sep_lists = ga_setup_selection.chunks(energies, 2)#splits up the list of results into the correct calculations
    for molen in sep_lists:#the rest is just the numbers from the orig 4 calculations correctly summed
        print molen
        x = (float(molen[1][1][1]) - float(molen[1][0][1]))
        y =  (float(molen[0][1][1]) - float(molen[0][0][1]))
        ea_ha = (float(molen[0][0][1]) - float(molen[1][0][1]))#this is electron affinity
        re_ha = x + y
        re_ev = re_ha *27.212
        ea_ev = ea_ha *27.212
        result.append((molen[0][0][0], re_ev, ea_ev))
    return result


def reorg_en1(all_mol_coords):
    '''runs initial reorganisation energy calculations'''
    com_files_neut = ga_setup_selection.com_generator_neutral_grndst(all_mol_coords)#the com files need to be made for both negative and neutral calcs
    com_files_neg = ga_setup_selection.com_generator_neg_mol_neg_geom(all_mol_coords)
    com_files = com_files_neut + com_files_neg#but they all are ran the same so the lists are combined
    job = group([run_gauss.s(com) for com in com_files])#run_gauss is a celery task, the .s is what allows the task to be passed to iridis
    run_job = job.apply_async()#this is what runs the job of all the com files, none are returned till they're all done
    result1 = run_job.get()#this returns the results
    return result1

def reorg_en2(opt_structs):
    '''second part of reorganisation energy calculations'''
    neg_structs = [(name,struct) for name,struct in opt_structs if "neg" in name]#this gets the mol name and its optimised structure 
    neut_structs = [(name,struct) for name,struct in opt_structs if "neg" not in name]
    com_files_neg = ga_setup_selection.com_generator_neg_mol_neut_charge(neg_structs)#creatng com files
    com_files_neut = ga_setup_selection.com_generator_neut_mol_neg_charge(neut_structs)
    com_files = com_files_neg + com_files_neut
    job = group([run_gauss2.s(com) for com in com_files])
    run_job = job.apply_async()
    energies = run_job.get()
    return energies

def reorg_en_calc(all_mol_coords):
    '''controls entire reorg en workflow'''
    result1 = reorg_en1(all_mol_coords)#do the first optimisations
    opt_structs = [(name,struct) for name,en,struct in result1 if (name,struct)]#get the structures to pass to the 2nd part
    first_energies = [(name,en) for name,en,struct in result1 if (name,en)]#get the energies for the sum
    second_energies = reorg_en2(opt_structs)#get the single point energies
    fe_sorted = sorted(first_energies, key=lambda tup: int(tup[0].split("_")[0]))#order first energies by name of mol
    se_sorted = sorted(second_energies, key=lambda tup: int(tup[0].split("_")[0]))#order second energies by name of mol
    energies = zip(fe_sorted,se_sorted)#join lists
    gen_result = reorg_en_calculator(energies)#calculate reorganisation energy
    return gen_result

def hl_diff_calc(all_mol_coords):#all the other calculations are much simpler and just involve one opt
    '''controls homo lumo difference workflow'''
    com_files = ga_setup_selection.com_generator_neutral_grndst(all_mol_coords)
    job =  group([run_gausshl.s(com) for com in com_files])
    run_job = job.apply_async()
    result1 = run_job.get()
    return result1

def lumo_calc(all_mol_coords):
    '''controls lumo energy workflow'''
    com_files = ga_setup_selection.com_generator_neutral_grndst(all_mol_coords)
    job = group([run_gaussl.s(com) for com in com_files])
    run_job = job.apply_async()
    result1 = run_job.get()
    return result1

def max_dipole_calc(all_mol_coords):
    '''controls the max dipole workflow'''
    com_files = ga_setup_selection.com_generator_neutral_grndst(all_mol_coords)
    job = group([run_gaussmd.s(com) for com in com_files])
    run_job = job.apply_async()
    result1 = run_job.get()
    return result1

