#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
author
------
  Novimir Antoniuk Pablant
    - npablant@pppl.gov
    - novimir.pablant@amicitas.com

description
-----------
  Launch a set of sfincs runs based on a runspec.dat file.
  Assume that these runs represent a scan in radius and radial electric field.
 
  This scan fullfills the need to perform radius and radial electric field
  scans that are specified externally.  In many ways this will serve the same
  purpose as sfincsScan 5, only specified differently.

  This scan will also put all of the sfincs runs into a single slurm batch
  job which can be launched as a single item.

  For the time being I am going to mimic the basic working of the sfincsScan
  scripts, though my hope is to eventually do a significant rewrite of
  the whole scrips set to make everything more generalized and easier to
  read/understand/maintain.

todo
----
  - Add a check of the input file to be sure that the species are in the
    expected order.
"""

import numpy as np
import logging
import os

from collections import OrderedDict

import mirscan_const
import mirscan_util

from mirscan_const import options


def get_job_name(param):
    return "sfincs"


def get_job_path(run):
    """
    Get a path for this job.

    This path is specific to mirscan_31 with the idea that I will use this 
    script to generate a SFINCS database.
    """
    #print('get_job_path')
    
    rN = mirscan_util.readVariable('rN_wish', 'float', run['inputfile'])
    
    nHats = mirscan_util.readVariable('nHats', 'float', run['inputfile'])
    dnHatdrHats = mirscan_util.readVariable('dnHatdrHats', 'float', run['inputfile'])
    THats = mirscan_util.readVariable('THats', 'float', run['inputfile'])
    dTHatdrHats = mirscan_util.readVariable('dTHatdrHats', 'float', run['inputfile'])
    
    Er = mirscan_util.readVariable('Er', 'float', run['inputfile'])
    
    #print('nHats: {}'.format(THats))
    #print('dnHatsdrHats: {}'.format(THats))    
    #print('THats: {}'.format(THats))
    #print('dTHatsdrHats: {}'.format(THats))
    
    path = ''
    
    rN_fmt = 'rN_{:0.3f}'    
    rN_path = rN_fmt.format(rN)
    path = os.path.join(path, rN_path)
    
    nT_fmt = 'n1_{:0.3f}_dndr1_{:0.3f}_t1_{:0.3f}_dtdr1_{:0.3f}_n2_{:0.3f}_dndr2_{:0.3f}_t2_{:0.3f}_dtdr2_{:0.3f}'
    nT_path = nT_fmt.format(
        nHats[0]
        ,dnHatdrHats[0]
        ,THats[0]
        ,dTHatdrHats[0]
        ,nHats[1]
        ,dnHatdrHats[1]
        ,THats[1]
        ,dTHatdrHats[1]
        )
    path = os.path.join(path, nT_path)
    
    Er_fmt = 'Er_{:0.3f}'
    Er_path = Er_fmt.format(Er)
    path = os.path.join(path, Er_path)
    
    logging.debug('Generated job path:')
    logging.debug(path)
    
    return path
    
    
def generate_combined_jobfiles_rn(run_list, jobfile_orig):
    # This routine is specific for SLURM.
    #
    # This is currently written so as to group all of the runs for a particular
    # radial position into a single job file.   This is particularly useful
    # for example when doing radial electric field scans.
    
    
    
    # First I need to split up the runs by radius.
    
    rN_list = list()
    for run in run_list:
        rN_list.append(mirscan_util.readVariable('rN_wish', 'float', run['inputfile']))
              
    rN_unique = np.unique(rN_list)
                       
    for ii_run, run in enumerate(run_list):
        run['job_index'] = np.flatnonzero(rN_unique == rN_list[ii_run])[0]


    jobfile_template = jobfile_orig.copy()                        
                            
    # Extract the line with the actual call to sfincs.
    num_matches = 0
    for ii_line, line in enumerate(jobfile_template):
        # skip empty and commented lines.
        if len(line.strip()) == 0:
            continue
        if line.strip()[0] in ('#', '!', '%'):
            continue
            
        if 'sfincs' in line:
            num_matches += 1
            ii_insert = ii_line
            sfincs_command = line
            del(jobfile_template[ii_line])

    if num_matches == 0:
        raise Exception('No call to sfincs found in job file.')
    if num_matches > 1:
        raise Exception('More than one call to sfincs found in job file.')
        
    # Create a set of jobfiles.    
    job_list = list()
    for rN in rN_unique:
        job_info = OrderedDict()
        job_info['rN'] = rN
        job_info['jobfile'] = jobfile_template.copy()
        job_info['path'] = os.path.join(options['path_base'], 'rN_{:0.3f}'.format(rN))
        job_info['filename'] = options['job_filename']+'_31'
        
        job_list.append(job_info)
                            
    
    # Insert all of the new sfincs commands into the combined jobfile at
    # the same location as the orginial command.
    for ii_run, run in enumerate(run_list):
        job_info = job_list[run['job_index']]
        jobfile = job_info['jobfile']
        
        path_full = os.path.join(options['path_base'], run['path'])            

        run_command = 'cd {} && {}\n'.format(path_full, sfincs_command)
        jobfile.insert(ii_insert, run_command)
        ii_insert += 1
        
        
    # Now write the files to disk.
    for job_info in job_list:
        filepath = os.path.join(job_info['path'], job_info['filename'])
        
        with open(filepath, "w") as ff:
            ff.writelines(job_info['jobfile'])
        logging.info('Wrote: {}'.format(filepath))
    
    
    return job_list

    
def generate_combined_jobfiles_all(run_list, jobfile):
    # This routine is specific for SLURM.
    #
    # This is currently written so that all of the jobs will be run in series.
    # It would be straight forward to modify this so that groups of jobs were
    # run in parallel.
    
    jobfile_all = jobfile.copy()
        
    # Extract the line with the actual call to sfincs.
    num_matches = 0
    for ii_line, line in enumerate(jobfile_all):
        # skip empty and commented lines.
        if len(line.strip()) == 0:
            continue
        if line.strip()[0] in ('#', '!', '%'):
            continue
            
        if 'sfincs' in line:
            num_matches += 1
            ii_insert = ii_line
            sfincs_command = line
            del(jobfile_all[ii_line])

    if num_matches == 0:
        raise Exception('No call to sfincs found in job file.')
    if num_matches > 1:
        raise Exception('More than one call to sfincs found in job file.')
    
    
    # Insert all of the new sfincs commands into the combined jobfile at
    # the same location as the orginial command.
    for ii_run, run in enumerate(run_list): 
        path_full = os.path.join(options['path_base'], run['path'])            

        run_command = 'cd {} && {}\n'.format(path_full, sfincs_command)
        jobfile_all.insert(ii_insert, run_command)
        ii_insert += 1
        
    jobfile_all_filename = options['job_filename']+'_31'
    jobfile_all_filepath = os.path.join(options['path_base'], jobfile_all_filename)
    with open(jobfile_all_filepath, "w") as ff:
        ff.writelines(jobfile_all)
    logging.info('Wrote: {}'.format(jobfile_all_filepath))
        
    
    output = list()
    output.append({
        'filepath':jobfile_all_filepath
        ,'filename':jobfile_all_filename
        })
    
    return output


def generate_jobarrayfile(run_list, jobfile):
    jobarrayfile_filename = options['jobarray_filename']
    jobarrayfile_filepath = os.path.join(options['path_base'], jobarrayfile_filename)
    
    output = list()
    output.append({
        'filepath':jobarrayfile_filepath
        ,'filename':jobarrayfile_filename
        })
    
    return output


def generate_jobarray_list(run_list):
    # This routine is specific for SLURM.
    #
    # Create an file than can be used to start an produce a slurm array.

    jobfile_out = []
    
    # Insert all of the new sfincs commands into the combined jobfile at
    # the same location as the orginial command.
    for ii_run, run in enumerate(run_list): 
        path_full = os.path.join(options['path_base'], run['path'])            
        jobfile_out.append(path_full+'\n')
        
    jobfile_out_filename = 'job_array.txt'
    jobfile_out_filepath = os.path.join(options['path_base'], jobfile_out_filename)
    with open(jobfile_out_filepath, "w") as ff:
        ff.writelines(jobfile_out)
    logging.info('Wrote: {}'.format(jobfile_out_filepath))
        
    
    output = list()
    output.append({
        'filepath':jobfile_out_filepath
        ,'filename':jobfile_out_filename
        })
    
    return output

                
def startscan(user_options=None):
    logging.debug('Entering mirscan_31.')
    
    mirscan_util.set_defaults_env()
    
    print(options)
    
    # Load the input file:
    with open(options['input_filename'], 'r') as ff:
        inputfile = ff.readlines()
        
    scan_type = mirscan_util.readScanVariable('scanType', 'int', inputfile)
    
    if not scan_type == 31:
        raise Exception(
            'scan_type is not equal to 31. Use generic mirscan command to launch'
            ' appropriate scan.')
    
    options['runspec_filename'] = mirscan_util.readScanVariable(
        "runSpecFile"
        ,"string"
        ,inputfile
        ,required=False
        ,stringValueCaseSensitive=True)
    
    if options['runspec_filename'] == None :
        options['runspec_filename'] = "runspec.dat"
    
    runspec_list = mirscan_util.readRunspec(options['runspec_filename'], verbose=True)
    mirscan_util.check_runspec_list(runspec_list, inputfile)
    
    
    # Create a more general run list to use in generating our runs.
    run_list = []
    for runspec in runspec_list:
        run = {}
        run['runspec'] = runspec
        run_list.append(run)
        
    # Read in the job.sfincsScan file:
    with open(options['job_filename'], 'r') as ff:
        jobfile = ff.readlines()
        
    # Read in the job_array.sfincsScan file:
    with open(options['jobarray_filename'], 'r') as ff:
        jobarrayfile = ff.readlines()
        
    # Loop over each entry found in the runspec.dat and generate new 
    # inputfile and jobfile lists. Also generate the paths for these runs. 
    for ii, run in enumerate(run_list):
        run['runspec']['name'] = get_job_name(runspec)
        
        logging.info("Beginning to handle job "+str(ii)+" of "+str(len(run_list))+": "+run['runspec']['name'])
                
        run['jobfile'] = mirscan_util.patch_jobfile(run['runspec'], jobfile)
        run['inputfile'] = mirscan_util.patch_inputfile(run['runspec'], inputfile) 
            
        run['path'] = get_job_path(run)

    
    # Create the base path if does not exist (this will usually be the current
    # working directory, so nothing will be done).
    mirscan_util.create_path(options['path_base'])
        
    # Create all of the input and job files in their final directories.
    # Currently I am not doing any check before overwriting.
    for ii, run in enumerate(run_list):
        
        path_full = os.path.join(options['path_base'], run['path'])
        mirscan_util.create_path(path_full)
        
        # Write the new inputfile.
        inputfile_filepath = os.path.join(path_full, options['input_filename'])
        with open(inputfile_filepath, "w") as ff:
            ff.writelines(run['inputfile'])
        logging.info('Wrote: {}'.format(inputfile_filepath))
            
        # Write the new jobfile.
        jobfile_filepath = os.path.join(path_full, options['job_filename'])
        with open(jobfile_filepath, "w") as ff:
            ff.writelines(run['jobfile'])
        logging.info('Wrote: {}'.format(jobfile_filepath))
    
        
    # Write the new jobarray file.
    
    # Setup the SLURM array command.
    # This is not a good place for this. At the very least this needs to be
    # in a separate function.
    #
    # WARNING: Number of steps per job current hardcoded.
    #          I need to find a way to automate this in some reasonable way.
    #
    #          I should have a few ways to control this.
    #            1. User option
    #            2. Specified in the input file as a !ss directive.
    #            3. Generated based on factorization of the total number
    #               runs.  This bit of code will return all factors.
    #               test = np.arange(n)+1
    #               test[np.where(n%(test) == 0)[0]]
    # 
    #          I'll need to replace a second line in the job_array file where
    #          number of runs per job is specified (RUNSPERJOB).
    runsperjob = 29
    numjobs = len(run_list)
    numjobarrays = int(np.ceil(numjobs/runsperjob))
    
    jobarray_options = {}
    jobarray_options['RUNSPERJOB'] = runsperjob
    jobarray_options['array'] = '1-{:0d}'.format(numjobarrays)
    jobarrayfile = mirscan_util.patch_jobfile(jobarray_options, jobarrayfile)
        
    jobarrayfile_filepath = os.path.join(options['path_base'], options['jobarray_filename'])
    with open(jobarrayfile_filepath, "w") as ff:
        ff.writelines(jobarrayfile)
    logging.info('Wrote: {}'.format(jobarrayfile_filepath))
            
    
    # Now generate the combined all job file.
    # job_runs = generate_combined_jobfiles_all(run_list, jobfile)
    
    # Now generate the combined rN job files.
    #job_runs = generate_combined_jobfiles_rn(run_list, jobfile)
     
    # Generate a files that can be used to create a slurm job array.
    job_paths = generate_jobarray_list(run_list)
    job_runs = generate_jobarrayfile(run_list, jobarrayfile)
    
    while True:
        proceed = input("Should I go ahead and launch these "+str(len(job_runs))+" jobs? [y/n] ")
        proceed = proceed.lower()
        if proceed[0] == "y" or proceed[0] == "n":
            break
        print("You must enter either y or n.")
    
    if proceed == "n":
        return None

    logging.info("launching jobs...")        
    
    mirscan_util.submit_jobs(job_runs)


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    #logging.basicConfig(level=logging.INFO)
    
    
    startscan()
    exit(0)
