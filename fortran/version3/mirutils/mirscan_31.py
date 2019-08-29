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
    
    
def startscan(user_options=None):
    logging.debug('Entering mirscan_31.')
    
    mirscan_util.set_defaults_env()
    
    print(options)
    
    # Load the input file:
    with open(options['input_filename'], 'r') as ff:
        inputfile_list = ff.readlines()
        
    scan_type = mirscan_util.readScanVariable('scanType', 'int', inputfile_list)
    
    if not scan_type == 31:
        raise Exception(
            'scan_type is not equal to 31. Use generic mirscan command to launch'
            ' appropriate scan.')
    
    options['runspec_filename'] = mirscan_util.readScanVariable(
        "runSpecFile"
        ,"string"
        ,inputfile_list
        ,required=False
        ,stringValueCaseSensitive=True)
    
    if options['runspec_filename'] == None :
        options['runspec_filename'] = "runspec.dat"
    
    runspec_list = mirscan_util.readRunspec(options['runspec_filename'], verbose=True)
    mirscan_util.check_runspec_list(runspec_list, inputfile_list)
    
    
    # Create a more general run list to use in generating our runs.
    run_list = []
    for runspec in runspec_list:
        run = {}
        run['runspec'] = runspec
        run_list.append(run)
        
    # Read in the job.sfincsScan file:
    with open(options['job_filename'], 'r') as ff:
        jobfile_list = ff.readlines()
        
    # Loop over each entry found in the runspec.dat and generate new 
    # inputfile and jobfile lists. Also generate the paths for these runs. 
    for ii, run in enumerate(run_list):
        run['runspec']['name'] = get_job_name(runspec)
        
        logging.info("Beginning to handle job "+str(ii)+" of "+str(len(run_list))+": "+run['runspec']['name'])
                
        run['jobfile'] = mirscan_util.patch_jobfile_list(run['runspec'], jobfile_list)
        run['inputfile'] = mirscan_util.patch_inputfile_list(run['runspec'], inputfile_list) 
            
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

            
    # Now generate the combined job file.
    # This is currently being written in-line and only for the SLURM manager.
    # I should move this into separate functions and also make sure that it is
    # general for all work managers.
    #
    # This is currently written so that all of the jobs will be run in series.
    # It would be straight forward to modify this so that groups of jobs were
    # run in parallel.
                
    jobfile_all = jobfile_list.copy()
        
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
        
            
        
        
        
        # Submit the sfincs job:
        #try:
        #    # We need to include .split(" ") to separate the command-line arguments into an array of strings.   
        #    submissionResult = subprocess.call(submitCommand.split(" "))
        #except:
        #    logging.error("ERROR! Unable to submit run "+jobName+" for some reason.")
        #    raise
        #else:
        #    if submissionResult==0:
        #        logging.info("No errors submitting job "+jobName)
        #    else:
        #        logging.error("Nonzero exit code returned when trying to submit job "+jobName)
        #
        #os.chdir("..")

        
        
        # dirNum = run_num - 1
        # If run directory already exists then skip this parameter.
        #
        # TODO:  Here is where I want to create a new directory structure.
        #while True:
        #    dirNum += 1
        #    job_name = str(dirNum)
        #    if dirNum < 10:
        #        job_name = "0" + job_name
        #    if not os.path.exists(job_name):
        #        break
        #os.mkdir(job_name)
        #os.chdir(job_name)
    
        #with open(options['job_filename'], "w") as ff:
        #    f.write(jobfile_new)
    
        #with open(options['input_filename'], "w") as ff:
        #    f.write(inputfile_new)        
        
    # I am going to separate the process of creating input files and launching jobs.
    # This question should be moved to the point where the jobs are actually launched.
    #while True:
    #    proceed = input("Should I go ahead and launch these "+str(num_runs_in_scan)+" jobs? [y/n] ").lower()
    #    if proceed[0] == "y" or proceed[0] == "n":
    #        break
    #    print("You must enter either y or n.")
    #
    #if proceed=="n":
    #    exit(0)
    #print("launching jobs...")
    
        
if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    
    startscan()
    exit(0)
