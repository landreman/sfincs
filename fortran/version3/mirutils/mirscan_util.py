# -*- coding: utf-8 -*-
"""
author
------
  Novimir Antoniuk Pablant
    - npablant@pppl.gov
    - novimir.pablant@amicitas.com

description
-----------
  A set of utilites to be used for the sfincsscan modules.

"""

import numpy as np
import logging
import os
import re
import subprocess

import mirscan_const
from mirscan_const import options

def create_path(path):
    """
    Check if the given path exits.  If not create the derectory and all 
    necessary parents.
    """
    if not os.path.exists(path):
        logging.info('Creating directory: {}'.format(path))
        os.makedirs(path)
        
        
def set_defaults_env():
    """
    Check the system environment for proper installation.

    todo:
      This routine is currently a mix of setting defaults, reading user options 
      and checking inputs.  These tasks should all be separated eventually.
    """
    
    options['path_base'] = os.getcwd()
    options['path_current'] = os.getcwd()
    
    if not os.path.isfile(options['input_filename']):
        logging.error(
            "Error! The file "+options['input_filename']+
            " must be present in the directory from which you call sfincsScan."
            )
        raise Exception('Input file not found.')
    
    try:
        options['sfincs_system'] = os.environ['SFINCS_SYSTEM']
    except:
        raise Exception("Error! Unable to read the SFINCS_SYSTEM environment variable. Make sure you have set it.")
    
    logging.info("I detect SFINCS_SYSTEM = "+options['sfincs_system'])
    
    # Define the Workload manager for each system.
    if options['sfincs_system'] == 'edison':
        options['job_manager'] = 'slurm'
    elif options['sfincs_system'] == 'cori':
        options['job_manager'] = 'slurm'
    elif options['sfincs_system'] == 'hydra':
        options['job_manager'] = 'loadleveler'
    elif options['sfincs_system'] == 'draco':
        options['job_manager'] = 'slurm'
    elif options['sfincs_system'] == 'portal':
        options['job_manager'] = 'slurm'
    elif options['sfincs_system'] == 'eddy':
        options['job_manager'] = 'slurm'  
    elif options['sfincs_system'] == 'laptop':
        options['job_manager'] = 'bash'
    elif options['sfincs_system'] == 'macports':
        options['job_manager'] = 'bash'
    elif options['sfincs_system'] == 'ubuntu16.04':
        options['job_manager'] = 'bash'
    else:
        logging.error("Error! SFINCS_SYSTEM="+options['sfincs_system']+" is not yet recognized by sfincsScan")
        logging.error("You will need to edit mirscan_util to specify a few things for this system.") 
        raise Exception("System defined in SFINCS_SYSTEM not recognized.")

    if not os.path.isfile(options['job_filename']):
        logging.error(
            "Error! A "+options['job_filename']+" file must be present in the directory"
            +" from which you call sfincsScan (even for systems with no queue)."
            )
        logging.error("Examples are available in /fortran/version3/utils/job.sfincsScan.xxx")
        raise Exception("Job file not found.")
    

def readScanVariable(
        varName
        ,intOrFloatOrString
        ,inputfile=None
        ,required=True
        ,stringValueCaseSensitive=False
        ):
    """
    This subroutine reads the special scan commands in the input.namelist that 
    are hidden from fortran.

    Todo:
      This routine is exactly the same as readVariable, except that it is looking
      lines that are preceded with '!ss'.  The readVariable routine should be
      generalized to handle this, and then this routine removed (or modifed to
      simply call readVariable with the right keywords.
    """
    if inputfile is None:
        raise NotImplementedError('Currently the input.namelist list must be given to this function.')

    if (intOrFloatOrString != 'int') and (intOrFloatOrString != 'float') and (intOrFloatOrString != 'string'):
        logging.error("intOrFloatOrString must be int, float, or string.")
        exit(1)

    originalVarName = varName
    varName = varName.lower()
    returnValue = None
    numValidLines = 0
    for line in inputfile:
        if not stringValueCaseSensitive:
            line2 = line.strip().lower()
        else:
            line2 = line.strip()

        # We need enough characters for the comment code, varName, =, and value:        
        if len(line2)<len(options['comment_code'])+3:
            continue

        if not line2[:len(options['comment_code'])]==options['comment_code']:
            continue

        line3 = line2[len(options['comment_code']):].strip()

        if len(line3) < len(varName)+2:
            continue

        if not line3[:len(varName)].lower()==varName:
            continue

        line4 = line3[len(varName):].strip()

        if not line4[0] =='=':
            continue

        line5 = line4[1:].strip();

        if intOrFloatOrString != 'string':
            # python does not recognize fortran's 1d+0 scientific notation
            line5 = line5.replace('d','e').replace('D','e')
        
        # Remove any comments:
        if '!' in line5:
            line5 = line5[:string.find(line5,'!')]
        line5 = line5.strip();

        if intOrFloatOrString=='int':
            try:
                returnValue = int(line5)
                numValidLines += 1
            except:
                logging.warning("Warning! I found a definition for the variable "+originalVarName+" in "+options['input_filename']+" but I was unable to parse the line to get an integer.")
                logging.warning("Here is the line in question:")
                logging.warning(line)
        elif intOrFloatOrString=='float':
            try:
                returnValue = float(line5)
                numValidLines += 1
            except:
                logging.warning("Warning! I found a definition for the variable "+originalVarName+" in "+options['input_filename']+" but I was unable to parse the line to get a float.")
                logging.warning("Here is the line in question:")
                logging.warning(line)
        elif intOrFloatOrString=='string':
            returnValue = line5
            numValidLines += 1

    if required and returnValue==None:
        logging.error("Error! Unable to find a valid setting for the scan variable "+originalVarName+" in "+options['input_filename']+".")
        logging.error("A definition should have the following form:")
        if intOrFloatOrString == 'int':
            logging.error(options['comment_code']+' '+originalVarName+' = 1')
        elif intOrFloatOrString == 'float':
            logging.error(options['comment_code']+' '+originalVarName+' = 1.5')
        elif intOrFloatOrString == 'string':
            logging.error(options['comment_code']+' '+originalVarName+' = nuPrime')
        exit(1)

    if numValidLines > 1:
        logging.warning("Warning! More than 1 valid definition was found for the variable "+originalVarName+". The last one will be used.")

    logging.info("Read "+originalVarName+" = "+str(returnValue))
    return returnValue


def parse_namelist_line(line, data_type=None):
    """
    Parse a singleline from a namelist file.  
    This will not work for multi-line entries.
    """
    
    if data_type is None: data_type = 'string'
        
    if data_type == 'int':
        convert_func = int
    elif data_type == 'float':
        convert_func = float
    elif data_type == 'string':
        convert_func = str
    else:
        raise Exception("data_type must be int, float, or string.")
    
    m = re.match(r'\s*(\w*)\s*=\s*([\d.eEdD,\- \t]*)', line)
    
    varname = m.group(1)
    
    values_string = m.group(2)
    values_string = values_string.lower().replace('d','e').strip()
    values_string_list = re.split(r'[\s,]+', values_string)
    try:
        values_list = [convert_func(xx) for xx in values_string_list]
    except:
        logging.warning('Variable {} found by conversion to {} failed.'.format(var_name, data_type))
        values_list = [None]
        
    return varname, values_list
    
    
def readVariable(var_name, data_type, inputfile=None, required=True):
    # This function reads normal fortran variables from the input.namelist file.
    # It is assumed that the input.namelist file has been loaded into the variable "inputFile".

    logging.debug('Searching for {} in inputfile.'.format(var_name))
        
    if (data_type != 'int') and (data_type != 'float') and (data_type != 'string'):
        raise Exception("data_type must be int, float, or string.")

    if data_type == 'int':
        convert_func = int
    elif data_type == 'float':
        convert_func = float
    elif data_type == 'string':
        convert_func = str
    else:
        raise Exception("data_type must be int, float, or string.")
    
    num_matches = 0
    for line in inputfile:
        m = re.match(r'\s*{}\s*=\s*([\d.eEdD,\- \t]*)'.format(var_name), line)
        if m:
            num_matches += 1
            values_string = m.group(1)
            values_string = values_string.lower().replace('d','e').strip()
            values_string_list = re.split(r'[\s,]+', values_string)
            try:
                values_list = [convert_func(xx) for xx in values_string_list]
            except:
                logging.warning('Variable {} found by conversion to {} failed.'.format(var_name, data_type))
                values_list = [None]
                    
    if num_matches == 0 and required:
        raise Exception(
            "Unable to find a valid setting for "+var_name+" in inputfile."
            )
    
    if num_matches > 1:
        logging.warning(
            "More than 1 valid definition was found for "+var_name+". The last one will be used."
            )
    
    if num_matches > 0:
        if len(values_list) == 1:
            output = values_list[0]
        else:
            output = values_list
        
        logging.debug("Read "+var_name+" = "+str(output)+" from inputfile.")
    else:
        output = None
        
    return output


def readDefault(varName, intOrFloatOrString, required=True):
    # This function reads the default value of fortran variables defined in globalVariables.F90.
    # If found it returns the last occurence of the variable, otherwise None.

    if (intOrFloatOrString != 'int') and (intOrFloatOrString != 'float') and (intOrFloatOrString != 'string'):
        logging.error("intOrFloatOrString must be int, float, or string.")
        exit(1)

    originalVarName = varName
    #varName = varName.lower()                                                                                                                      
    returnValue = None
    numValidLines = 0

    try: 
        working_dir = os.getcwd() ##Store current working directory
        os.chdir(os.path.dirname(os.path.abspath(__file__))) ##Go to directory of this file
        defaultVariablesFile = open('../' + defaultVariablesFilename, 'r') ##Open file
        os.chdir(working_dir) ##Go back to working directory
    except:
        logging.error("Error! Unable to open "+defaultVariablesFilename+".")
        if required:
            raise
        else:
            return returnValue

    for line in defaultVariablesFile:

        #line3 = line.strip().lower()                                                                                                               
        line3 = line.strip()
        if len(line3)<1:
            continue

        if line3[0]=='!':
            continue

        if len(line3) < len(varName)+2:
            continue

        #if not line3[:len(varName)].lower()==varName.lower():
        #    continue

        begin_index = line3.lower().find(varName.lower())
        
        if begin_index == -1: #Cannot find varName on this line
            continue

        if begin_index != 0 and line3[begin_index-1] != ' ': #If character before varName is not a blank space, this is the wrong variable  
            continue

        line3 = line3.replace(' ', '')

        start_index = line3.lower().find(varName.lower())
                
        line4 = line3[start_index:].strip()

        
        if len(line4) < len(varName)+2:
            continue 

        if not line4[len(varName)] =='=':
            continue

        line5 = line4[len(varName)+1:].strip();
        line5 = line5.split(',')[0] ##Needed if several variables are defined on the same line  

        if intOrFloatOrString != 'string':
            # python does not recognize fortran's 1d+0 scientific notation                                                                          
            line5 = line5.replace('d','e').replace('D','e')
 
        # Remove any comments:                                                                                                                      
        if '!' in line5:
            line5 = line5[:string.find(line5,'!')]

           

        if intOrFloatOrString=='int':
            try:
                returnValue = int(line5)
                numValidLines += 1
            except:
                logging.warning(
                    "Warning! I found a definition for the variable "
                    +originalVarName+" in "+defaultVariablesFilename
                    +" but I was unable to parse the line to get an integer."
                    )
                logging.warning("Here is the line in question:")
                logging.warning(line)
        elif intOrFloatOrString=='float':
            try:
                returnValue = float(line5)
                numValidLines += 1
            except:
                logging.warning(
                    "Warning! I found a definition for the variable "
                    +originalVarName+" in "+defaultVariablesFilename
                    +" but I was unable to parse the line to get a float."
                    )
                logging.warning("Here is the line in question:")
                logging.warning(line)
        elif intOrFloatOrString=='string':
            returnValue = line5
            numValidLines += 1

    if required and returnValue==None:
        logging.error(
            "Error! Unable to find a valid setting for the variable "
            +originalVarName+" in "+defaultVariablesFilename+"."
            )
        exit(1)

    if numValidLines > 1:
        logging.warning(
            "Warning! More than 1 valid definition was found for the variable "
            +originalVarName+". The last one will be used."
            )
        
    logging.info("Read "+originalVarName+" = "+str(returnValue))
    return returnValue


def numberatend(string):
    """
    Extract the number at the end of a variable name.

    ToDo:
      - Rename to: extract_var_number_spec
      - Generalize using regex. Currently this implemenation will only work 
        up to 9 variables.
    """
    thenumb=0
    if string[-2]=='_':
        try:
            thenumb=int(string[-1])
        except:
            thenumb=0
    return thenumb


def readRunspec(runspec_filename, verbose=False):
    """
    Read runspect.dat formatted files for sfincsscan.

    ToDo:
      I should split parsing of the file from reading of the file.
      That would make it easier to deal with strings.
    """
    
    if not os.path.isfile(runspec_filename):
        raise Exception(
            "Error! The file "+runspec_filename
            +" must be present in the directory from which you call sfincsScan"
            +", or the full path must be given."
            )

    output = []
    
    # Load the runspec file:
    with open(runspec_filename, 'r') as ff:
        runspec_string = ff.readlines()
    
    # First read the file header.
    for ii_line, line in enumerate(runspec_string):
        if (
            runspec_string[ii_line][0] == '!' 
            or runspec_string[ii_line][0] == '%' 
            or runspec_string[ii_line][0] == '#'
            ):
            
            continue
        
        # Here I assume that there is only a single comment character.
        runparamnames = runspec_string[ii_line-1][1:].split()
        noElems = len(runparamnames)
        break
    
    # Now read each individual rine in the runspec file.
    for ii_line, line in enumerate(runspec_string):
        if (
            runspec_string[ii_line][0] == '!' 
            or runspec_string[ii_line][0] == '%' 
            or runspec_string[ii_line][0] == '#'
            ):
            
            continue

        # Skip empty lines.
        if len(line) <= 1:
            continue
        
        thesplittedstring = line.split()
        if not len(thesplittedstring) == noElems:
            raise Exception(
                "Error! In the runspec file, the number of variable"
                " names does not match the number of variables!"
                )

               
        entry = {} 
        for ii_par in range(noElems):
            name = runparamnames[ii_par]
            thestring=thesplittedstring[ii_par]
            thestring=thestring.replace('d','e').replace('D','e')
            if ('.' in thestring) or ('e' in thestring):
                value = float(thestring)
            else:
                value = int(thestring)
            
            entry[name] = value 
    
        output.append(entry)

        
    if verbose:
        # Print the data from runspec.dat to the screen
        print('runspec file read from:')
        print('  {}'.format(os.path.abspath(runspec_filename)))
        print(' '+' '.join(['{:>12s}'.format(v) for v in output[0].keys()]))
        for entry in output:
            oneline = ""
            for key, value in entry.items():
                theformat = "{:12.3e}"
                oneline = oneline+' '+theformat.format(value)
            print(oneline)
        
    return output


def check_runspec_list(runspec_list, inputfile):
    
    jm_strings = [
        "name"
        ,"node"
        ,"tasks_per_node"
        ,"wall_clock_limit"
        ,"nodes"
        ,"ntasks_per_node"
        ,"time"
        ]
    
    ##Check that variables are defined in SFINCS input file
    for paramname in runspec_list[0].keys():
        if len(str(paramname).split("_")) != 1 :
            try: 
                # Check if the last part of the variable name can be transformed into an integer
                # If this fails an exception is raised, then do nothing to paramname.
                dummy = int(str(paramname).split("_")[-1])
                paramname = '_'.join(str(paramname).split("_")[:-1])
            except :
                pass
            
        test_read_name = readVariable(paramname, "string", inputfile, required=False)
        if test_read_name == None:
            if not (paramname in jm_strings):
                logging.error(
                    "The variable " + paramname + 
                    " you wish to scan must be explicitly assigned some value"
                    " in the appropriate namelist, even though the value will"
                    " be ignored in the scan."
                    )
                raise Exception('runlist param must be defined in the input namelist.')
        

def patch_jobfile(params, jobfile_orig):
    
    # WARNING! This routine is currently broken.
    
    jobfile_new = jobfile_orig
 
    if options['job_manager'] == 'slurm':
        jobfile_new = patch_jobfile_slurm(params, jobfile_orig)
    elif options['job_manager'] == 'slurm':
        jobfile_new = patch_jobfile_loadleveler(params, jobfile_orig)
    elif options['job_manager'] == 'bash':
        jobfile_new = jobfile_orig.copy()
        
    return jobfile_new


def patch_jobfile_slurm(params, jobfile_orig):
    """
    todo
    ----
      This should accept either a list or a string.
    """
        
    jobfile_new = jobfile_orig.copy()
    
    # This next function is defined separately for each system in sfincsScan
    # mirscan_util.nameJobFile(thisJobFile, jobName)
    
    # Still need set the job name!
       
    for key, value in params.items():

        if key == "name":
            for ii, line in enumerate(jobfile_new):
                if "#SBATCH -J" in line:
                    jobfile_new[ii] = "#SBATCH -J "+str(value)+"\n"
                    
        if key == "nodes":
            for ii, line in enumerate(jobfile_new):
                if "#SBATCH --nodes" in line:
                    jobfile_new[ii] = "#SBATCH --nodes="+str(value)+"\n"

        if key == "ntasks-per-node":
            for ii, line in enumerate(jobfile_new):
                if "#SBATCH --ntasks-per-node" in line:
                    jobfile_new[ii] = "#SBATCH --ntasks-per-node="+str(value)+"\n"
                    
        if key == "array":
            for ii, line in enumerate(jobfile_new):
                if "#SBATCH --array" in line:
                    jobfile_new[ii] = "#SBATCH --array="+str(value)+"\n"
                    
        if key == "RUNSPERJOB":
            for ii, line in enumerate(jobfile_new):
                if "RUNSPERJOB=" in line:
                    jobfile_new[ii] = "RUNSPERJOB={:0d}\n".format(value)
                    
        if key == "time":
            hours = value//60
            minutes = value%60
            hour_str = '{:02d}'.format(hours)
            minute_str = '{:02d}'.format(minutes)
            for ii, line in enumerate(jobfile_new):
                if "#SBATCH --time" in line:
                    jobfile_new[ii] = "#SBATCH --time="+hour_str+":"+minute_str+":00\n"
 
    return jobfile_new


def patch_jobfile_loadleveler(params, jobfile_orig):
    """
    todo
    ----
      This should accept either a list or a string.
    """
        
    jobfile_new = jobfile_orig.copy()
        
    # Still need set the job name!
    
    for key, value in params.items():
        if key == "name":
            for ii, line in enumerate(jobfile_new):
                if "# @ job_name =" in line:
                    jobfile_new[ii] = "# @ job_name = "+str(value)+"\n"
                    
        if key == "node":
            for ii, line in enumerate(jobfile_new):
                if "# @ node =" in line:
                    jobfile_new[ii] = "# @ node = "+str(value)+"\n"

        if key == "tasks-per-node":
            for ii, line in enumerate(jobfile_new):
                if "# @ tasks_per_node =" in line:
                    jobfile_new[ii] = "# @ tasks_per_node = "+str(value)+"\n"

        if key == "ConsumableCpus":
            for ii, line in enumerate(jobfile_new):
                if "# @ resources = ConsumableCpus" in line:
                    jobfile_new[ii] = "# @ resources = ConsumableCpus("+str(value)+")\n"
                    
        if key == "wall_clock_limit":
            hours = value//60
            minutes = value%60
            hour_str = '{:02d}'.format(hours)
            minute_str = '{:02d}'.format(minutes)
            for ii, line in enumerate(jobfile_new):
                if "# @ wall_clock_limit =" in line:
                    jobfile_new[ii] = "# @ wall_clock_limit = "+hour_str+":"+minute_str+":00\n"

    return jobfile_new


def patch_inputfile(params, inputfile_orig):
    """
    Patch an inputfile list with a set of new parameters.

    ToDo:
      - This should accept either a list or a string.
    """
    inputfile_new = inputfile_orig.copy()
    
    for ii_line, line in enumerate(inputfile_new):
        if line[-22:] == " ! Set by sfincsScan.\n":
            # Remove this old scan note from the file. Only indicate the present scan variables.
            inputfile_new[ii_line] = line[:-22]+"\n"
            
        for key, value in params.items():
            postfix = numberatend(key)
            if postfix == 0:
                varname = key
            else:
                varname = key[:-2]
                
            if namelistLineContains(line, varname):
                varname_orig, values = parse_namelist_line(line)
            
                if postfix > len(values):
                    raise Exception('More values in runspec than in inputfile for {}'.format(varname_orig))
            
                values[postfix-1] = value
                values_string = ' '.join([str(xx) for xx in values])
                
                # It is important that we modify line in place, as well as
                # update the new inputfile list.
                line =  "  {} = {} ! Set by sfincsScan.\n".format(varname_orig, values_string)
                inputfile_new[ii_line] = line

    return inputfile_new


def get_submit_command():
    """
    Return the appropriate command for the specified job manager.
    """
    
    if options['job_manager'] == 'slurm':
        command = 'sbatch'
    elif options['job_manager'] == 'loadleveler':
        command = 'llsubmit'
    elif options['job_manager'] == 'bash':
        command = 'bash'
        
    return command
    

def submit_jobs(job_list):
    """
    Submit the given jobs to the workload manager.

    For now this is only written for SLURM, but should eventually work for the 
    other options as well.
    """
    
    submit_command = get_submit_command()

    for job in job_list:
        filepath = job['filepath']
        
        command = [submit_command, filepath]
        command_string = ' '.join(command)
            
        try:
            command_string = command_string
            result = subprocess.run(command_string, shell=True)
            
            if result.returncode == 0:
                logging.info('Successfuly submitted job: {}'.format(command_string))
            else:
                logging.info('Job submission had errors: {}'.format(command_string))                
        except:
            logging.exception("ERROR! Unable to submit run "+command_string+" for some reason.")   


def uniq(seq): 
   checked = []
   for e in seq:
       if e not in checked:
           checked.append(e)
   return checked


def logspace(min,max,nn):
    if nn < 1:
        return []
    elif nn==1:
        return [min]

    if min <= 0:
        log.error("Error in logspace! min must be positive.")
        exit(1)
    if max <= 0:
        log.error("Error in logspace! max must be positive.")
        exit(1)
    return [math.exp(x/(nn-1.0)*(math.log(max)-math.log(min))+math.log(min)) for x in range(nn)]


def linspace(min,max,nn):
    if nn < 1:
        return []
    elif nn==1:
        return [min]
    return [x/(nn-1.0)*(max-min)+min for x in range(nn)]


def logspace_int(min,max,nn):
    return uniq(map(int,map(round,logspace(min,max,nn))))


def logspace_odd(min,max,nn):
    temp = map(int,logspace(min,max,nn))
    temp2 = []
    for x in temp:
        if (x % 2 == 0):
            temp2.append(x+1)
        else:
            temp2.append(x)
    return uniq(temp2)


def namelistLineContains(line,varName):
    line2 = line.strip().lower()
    varName = varName.lower()
    # We need enough characters for the varName, =, and value: 
    if len(line2)<len(varName)+2:
        return False

    if line2[0]=="!":
        return False

    nextChar = line2[len(varName)]
    if line2[:len(varName)]==varName and (nextChar==" " or nextChar=="="):
        return True
    else:
        return False

    
def namelistLineContainsSS(line,varName):
    # Same as namelistLineContains, but looking for !ss directives.
    line2 = line.strip().lower()
    varName = varName.lower()
    if len(line2)<len(options['comment_code']):
        return False

    if line2[:len(options['comment_code'])] != options['comment_code']:
        return False

    # If we got this far, the line must begin with !ss, so strip this part out.
    line2 = line2[len(options['comment_code']):].strip()

    # We need enough characters for the varName, =, and value: 
    if len(line2)<len(varName)+2:
        return False

    if line2[0]=="!":
        return False

    nextChar = line2[len(varName)]
    if line2[:len(varName)]==varName and (nextChar==" " or nextChar=="="):
        return True
    else:
        return False
