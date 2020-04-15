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
  Analyse the output from mirscan_31 and determine the ambipolar Er.
"""

import numpy as np
import collections
import logging
import os
import glob
import h5py

from collections import OrderedDict

from scipy.interpolate import PchipInterpolator
from scipy.interpolate import interp1d
from scipy.optimize import brentq

import mirhdf5

import mirscan_const
import mirscan_util

from mirscan_const import options

import mirscan_31


def run(path=None):
    """
    Start analysis for the given path.

    For now this just assumes that the given directory contains a parameter
    scan, with each step containing an er scan.  In the future we can determine
    the scan type more intelegently.
    """
    output = analyse_ambipolar_param_scan(path)
    return output
    
    
def _open_hdf5(file_in):
    if isinstance(file_in, str):
        try:
            obj_file = h5py.File(file_in, 'r')
        except:
            raise Exception("Unable to open "+file_in+" even though this file exists.")
    elif isinstance(file_in, h5py.File):
        obj_file = file_in
    else:
        raise Exception('Input must be a string or an h5py File object.')
    
    return obj_file
    

def check_run_success(file_in):
    obj_file = _open_hdf5(file_in)
    
    try:
        finished = obj_file["finished"]
    except KeyError:
        raise Exception("Run in "+obj_file.filename+" does not appear to have finished.")
    
    try:
        # Try reading a field that should definitely be present in the output 
        # file for any run that completed.
        dummy = obj_file["FSABFlow"][...]
    except:
        raise Exception("Unable to read "+obj_file.filename+" even though this file exists.")


    # Check for nonlinear convergence when Phi1 is included.
    Check_integerToRepresentTrue = (obj_file["integerToRepresentTrue"][...])
    Check_includePhi1 = (obj_file["includePhi1"][...] == Check_integerToRepresentTrue)
    if Check_includePhi1 :
        try:
            if (obj_file["didNonlinearCalculationConverge"][...] != Check_integerToRepresentTrue):
                raise Exception("The nonlinear solver in "+obj_file.filename+" does not appear to have converged.")
        except KeyError:
            raise Exception("The nonlinear solver in "+obj_file.filename+" does not appear to have finished.")
    
    return True
   

def _output_append_value(dict_in, key=None, value=None, label=None, species=None):
    if species is not None:
        key += '_{}'.format(species)
        
    if not key in dict_in['value']:
        dict_in['value'][key] = []
    if not key in dict_in['label']:
        dict_in['label'][key] = label
        
    if isinstance(value, (collections.Sequence, np.ndarray)):
        for vv in value:
            dict_in['value'][key].append(vv)
    else:
        dict_in['value'][key].append(value)
            

def extract_results(path):
    """
    Extract results from a single sfincs run.

    This function will create a simplified dictionary of results from
    the complete sfincs output.

    Programming Notes
    -----------------
      The use of _output_append_value here is not needed since this function
      only deals with a single run.  I'm going to keep using it though just
      for consistency.
    """
                
    output = OrderedDict()
    output['value'] = OrderedDict()
    output['label'] = OrderedDict()
    output['param'] = OrderedDict()
    output['info'] = OrderedDict()

    output['info']['output_type'] = 'run'
    output['info']['path'] = path
    
    filepath_o = os.path.join(path, options['output_filename'])
    filepath_i = os.path.join(path, options['input_filename'])

    output['info']['output_file'] = filepath_o  
    output['info']['input_file'] = filepath_i

    obj_file = _open_hdf5(filepath_o)
    try:
        check_run_success(obj_file)
    except:
        raise Exception('Run '+obj_file.filename+' could not be successfuly read.')
    

        
        
    # Extract some information from the input file.
    # When exctacting data for a scan, this only needs to be done once.
    # If I want to speed things up a little I could make this optional.
    
    with open(filepath_i, 'r') as ff:
        inputfile = ff.readlines()
        
    inputRadialCoordinate = mirscan_util.readVariable(
        "inputRadialCoordinate"
        ,"int"
        ,inputfile
        ,required=False)
    if inputRadialCoordinate == None:
        inputRadialCoordinate = mirscan_util.readDefault("inputRadialCoordinate", "int", inputfile)
    output['param']['inputRadialCoordinate'] = inputRadialCoordinate

    if inputRadialCoordinate==0:
        radiusName = "psiHat"
    elif inputRadialCoordinate==1:
        radiusName = "psiN"
    elif inputRadialCoordinate==2:
        radiusName = "rHat"
    elif inputRadialCoordinate==3:
        radiusName = "rN"
    else:
        raise Exception("Error! Invalid inputRadialCoordinate.")

    output['param']['radiusName'] = radiusName
    
    radius_wish = mirscan_util.readVariable(
        radiusName+"_wish"
        ,"float"
        ,inputfile
        ,required=False)    
    if radius_wish == None:
        radius_wish = mirscan_util.readDefault(radiusName+"_wish","float")
        
    key = radiusName+'_wish'
    label = 'Radius Wish'
    value = radius_wish
    _output_append_value(output, key, value, label)

    
    # The expression [()] converts from an h5py dataset to a numpy ndarray:
    integerToRepresentTrue = (obj_file["integerToRepresentTrue"][()])
    inputRadialCoordinateForGradients = obj_file["inputRadialCoordinateForGradients"][()]
    RHSMode = obj_file["RHSMode"][()]
    Nspecies = obj_file["Nspecies"][()]
    Zs = obj_file["Zs"][()]
    includePhi1 = (obj_file["includePhi1"][()] == integerToRepresentTrue)

    output['param']['inputRadialCoordinateForGradients'] = inputRadialCoordinateForGradients
    output['param']['RHSMode'] = RHSMode
    output['param']['Nspecies'] = Nspecies
    output['param']['Zs'] = Zs
    output['param']['includePhi1'] = includePhi1

    
    if RHSMode != 1 and RHSMode != 2 and RHSMode != 3:
        raise Exception("Error! sfincsScanPlot is not yet set up for RHSMode = "+str(RHSMode))              
    
    key = radiusName
    label = 'Radius'
    value = obj_file[key][()]
    _output_append_value(output, key, value, label)
  

    if inputRadialCoordinateForGradients==0:
        ErName = "dPhiHatdpsiHat"
    elif inputRadialCoordinateForGradients==1:
        ErName = "dPhiHatdpsiN"
    elif inputRadialCoordinateForGradients==2:
        ErName = "dPhiHatdrHat"
    elif inputRadialCoordinateForGradients==3:
        ErName = "dPhiHatdrN"
    elif inputRadialCoordinateForGradients==4:
        ErName = "Er"
    else:
        raise Exception("Error! Invalid inputRadialCoordinateForGradients.")
    
    output['param']['ErName'] = ErName
    
    key = ErName
    label = ErName
    value = obj_file[ErName][()]
    _output_append_value(output, key, value, label)

    if RHSMode > 1:
        raise NotImplementedError('Handling of RHSMode > 1 not implemented.')
        #results = []            
        #transportMatrix = obj_file["transportMatrix"][()]
        #for icol in range(transportMatrix.shape[1]):
        #    for irow in range(transportMatrix.shape[0]):
        #        label = "L"+str(irow+1)+str(icol+1)
        #        value = transportMatrix[irow,icol]
        #output.append(results)
        #radialCurrents.append(0)
    else:
        # Build the output for RHSMode = 1

        iteration = -1
        for ii_spec in range(Nspecies):
            
            key = 'nHats'
            label = 'nHats (species {})'.format(ii_spec)
            value = obj_file[key][()][ii_spec]
            _output_append_value(output, key, value, label, ii_spec+1)
            
            key = 'THats'
            label = 'THats (species {})'.format(ii_spec)
            value = obj_file[key][()][ii_spec]
            _output_append_value(output, key, value, label, ii_spec+1)                
            
            key = 'dnHatdrHat'
            label = 'dnHatdrHat (species {})'.format(ii_spec)
            value = obj_file[key][()][ii_spec]
            _output_append_value(output, key, value, label, ii_spec+1)
            
            key = 'dTHatdrHat'
            label = 'dTHatdrHat (species {})'.format(ii_spec)
            value = obj_file[key][()][ii_spec]
            _output_append_value(output, key, value, label, ii_spec+1)

            
            key = 'FSABFlow'
            label = 'FSABFlow (species {})'.format(ii_spec)
            value = obj_file[key][()][ii_spec,iteration]
            _output_append_value(output, key, value, label, ii_spec+1)                
                                                        
            key = 'particleFlux_vm_rHat'
            label = 'particleFlux vm_rHat (species {})'.format(ii_spec+1)
            value = obj_file[key][()][ii_spec,iteration]
            _output_append_value(output, key, value, label, ii_spec+1)
            
            key = 'heatFlux_vm_rHat'
            label = 'heatFlux vm_rHat (species {})'.format(ii_spec+1)
            value = obj_file[key][()][ii_spec,iteration]
            _output_append_value(output, key, value, label, ii_spec+1)
                
            if includePhi1:
                
                key = 'particleFlux_vd_rHat'
                label = 'particleFlux vd_rHat (species {})'.format(ii_spec+1)
                value = obj_file[key][()][ii_spec,iteration]
                _output_append_value(output, key, value, label, ii_spec+1)
                
                key = 'heatFlux_vd_rHat'
                label = 'heatFlux vd_rHat (species {})'.format(ii_spec+1)
                value = obj_file[key][()][ii_spec,iteration]
                _output_append_value(output, key, value, label, ii_spec+1)
                
                key = 'heatFlux_withoutPhi1_rHat'
                label = 'heatFlux withoutPhi1_rHat (species {})'.format(ii_spec+1)
                value = obj_file[key][()][ii_spec,iteration]
                _output_append_value(output, key, value, label, ii_spec+1)

 
        key = 'FSABjHatOverRootFSAB2'
        label = 'FSABjHatOverRootFSAB2'
        value = obj_file[key][()][iteration]
        _output_append_value(output, key, value, label)   
        
        key = 'FSABjHat'
        label = 'FSABjHat'
        value = obj_file[key][()][iteration]
        _output_append_value(output, key, value, label)                
        
        key = 'radialCurrent'
        label = 'Radial Current'
        if includePhi1:
            value = sum(Zs*obj_file['particleFlux_vd_rHat'][()][:,iteration])
        else:
            value = sum(Zs*obj_file['particleFlux_vm_rHat'][()][:,iteration])
        _output_append_value(output, key, value, label)

        
        # These next two values are apparently needed when NSpecies = 1
        # Since I don't know what they are, I am going to keep them
        # commented out for now.
               
        #key = "sources"
        #label = "Source 1"
        #value = obj_file["sources"][()][0,0,0]
        #_output_append_value(output, key, value, label, 1)
                    
        #key = "sources"
        #label = "Source 2"
        #value = obj_file["sources"][()][1,0,0]
        #_output_append_value(output, key, value, label, 2) 

    logging.debug("Successfully read run "+obj_file.filename)

    return output

    
def extract_results_scan(path):
    
    # This routines currently assumes that in the given path are a set of
    # directories that contain the Er scan runs.

    search_string = os.path.join(path,'**/'+options['output_filename'])
    logging.debug('Searching for runs in: {}'.format(search_string))
    filepaths_output = glob.glob(search_string, recursive=True)
    filepaths_output.sort()
    paths_scan = [os.path.dirname(xx) for xx in filepaths_output]
            
    output = OrderedDict()
    output['value'] = OrderedDict()
    output['label'] = OrderedDict()
    output['param'] = OrderedDict()
    output['info'] = OrderedDict()

    output['info']['output_type'] = 'scan'
    output['info']['path'] = path

    first = True
    for ii, path in enumerate(paths_scan):
        try:
            output_run = extract_results(path)
        except:
            logging.error('Could not load results in path {}'.format(path))
            continue

        # The params should be the same for all file in the scan.
        if first:
            first = False
            for key in output_run['param']:
                output['param'][key] = output_run['param'][key]
                     
        for key in output_run['param']:
            if np.any(output['param'][key] != output_run['param'][key]):
                raise Exception('Parameter {} is not conistent between runs in scan.'.format(key))

        for key in output_run['value']:
            value = output_run['value'][key]
            label = output_run['label'][key]
            _output_append_value(output, key, value, label)

    return output


def analyse_ambipolar_param_scan(path=None):
    """
    Loop through each radii from mirscan_31 and do an analysis.
    """

    if path is None: path = ''

    if not os.path.exists(path):
        raise Exception('Given path does not exist.')
    
    # This routines currently assumes the specific mirscan_31 directory structure. 

    # Find a list of all completed runs.
    search_string = os.path.join(path, '**/'+options['output_filename'])
    paths_output = glob.glob(search_string, recursive=True)
    paths_output.sort()
    paths_er = [os.path.dirname(xx) for xx in paths_output]
        
    # The check for completed sfincs runs happens inside of the extraction
    # routines. For now we just assume that all runs were successful.
    
    # Each of these paths should contain an Er scan.
    paths_scan = [os.path.dirname(xx) for xx in paths_er]
    paths_scan = np.unique(paths_scan)
         
    output = OrderedDict()
    output['value'] = OrderedDict()
    output['label'] = OrderedDict()
    output['param'] = OrderedDict()
    output['info'] = OrderedDict()
    
    output['info']['output_type'] = 'ambipolar'
    output['info']['path'] = path

    if len(paths_scan) == 0:
        raise Exception('No ambipolar Er scans found in directory.')
    
    num_success = 0
    for ii, path_run in enumerate(paths_scan):
        try:
            output_scan = analyse_ambipolar_scan(path_run)
        except:
            logging.exception('Exception determining ambipolar solution for {}'.format(path_run))
            continue

        # The params should be the same for all file in the scan.
        if ii == 0:
            for key in output_scan['param']:
                output['param'][key] = output_scan['param'][key]
                     
        for key in output_scan['param']:
            if output['param'][key] != output_scan['param'][key]:
                raise Exception('Parameter {} is not conistent between runs in scan.'.format(key))

        for key in output_scan['value']:
            value = output_scan['value'][key]
            label = output_scan['label'][key]
            _output_append_value(output, key, value, label)

        num_success += 1

    if num_success == 0:
        raise Exception('No Er scan analysis were successful.')
    
    filename = 'ambipolarSolutions.h5'
    filepath = os.path.join(path, filename)
    mirhdf5.dictToHdf5(output, filepath)
    logging.info('Wrote ambipolar parameter scan solutions to: {}'.format(filepath))    
    
    return output


def load_ambipolar_scan(path):

    # First check if the scan results are in the current folder.
    filepath = os.path.join(path, 'ambipolarSolutions.h5')

    # Now try to find if the results are in a nested folder.
    if not os.path.exists(filepath):
        search_string = os.path.join(path, '**/ambipolarSolutions.h5')
        filepath_list = glob.glob(search_string, recursive=True)
        filepath_list.sort()
        if len(filepath_list) > 1:
            raise Exception('More than one ambipolarSolutions.h5 found.')
        
        if len(filepath_list) == 0:
            filepath = None
        else:
            filepath = filepath_list[0]

    if filepath is not None:
        output = mirhdf5.hdf5ToDict(filepath)
        logging.info('Read ambipolar solutions from: {}'.format(filepath))
    else:
        output = analyse_ambipolar_scan(path)

    

    return output


def analyse_ambipolar_scan(path):

    output = extract_results_scan(path)

    if output['param']['RHSMode'] != 1:
        raise Exception('Ambipolar solution can only be done with RHSMode = 1.')
    
    ErName = output['param']['ErName']
    conversionFactorToEr = output['value']['Er'][0]/ output['value'][ErName][0]
    
    # Sort by Er:
    sort_index = np.argsort(output['value'][ErName])
    for key in output['value']:
        output['value'][key] = np.array(output['value'][key])[sort_index]
        
    logging.debug("Here comes radialCurrent")
    logging.debug(output['value']['radialCurrent'])
    
    numQuantities = len(output['value'])
    
    # ***************************************************
    # Solve for the ambipolar E_r
    # ***************************************************
    
    # The PchipInterpolator routine seems to return an error when there are 
    # only 2 elements to the input array.
    if len(output['value'][ErName]) < 3:
        whichInterpolator = 1
    else:
        whichInterpolator = 0
    
    ambipolaritySolveSucceeded = False
 
    if max(output['value']['radialCurrent']) > 0 and min(output['value']['radialCurrent']) < 0:
        try:
            # We have runs that straddle j_psi = 0, so we can find at least 1 ambipolar solution.
            logging.info("Solving for E_r.")
            
            # Note: For interpolation it does not work well to use scipy's "interp1d" function, 
            #       since this tends to give large overshoots for the typical data we get in an 
            #       Er scan (with sharp structure near Er=0).
            if whichInterpolator==0:
                interpolator = PchipInterpolator(output['value'][ErName], output['value']['radialCurrent'])
            else:
                interpolator = interp1d(output['value'][ErName], output['value']['radialCurrent'])
             
            # In the next line, 500 is arbitrary, and could be replaced by another large integer.
            Er_fine = np.linspace(min(output['value'][ErName]), max(output['value'][ErName]), num=500)
            radialCurrent_fine = interpolator(Er_fine)
            
            # Approximately find points where the sign of radial current flips:
            positiveCurrent = (radialCurrent_fine > 0)
            signFlips = (positiveCurrent[:-1] != positiveCurrent[1:])
            numRoots = sum([1 for xx in signFlips if xx])
            logging.info("Number of roots found for E_r: {}".format(numRoots))
            
            roots = []
            #for i in range(numRoots):
            for index,value in enumerate(signFlips):
                if value:
                    roots.append(brentq(interpolator,Er_fine[index],Er_fine[index+1]))
            # Convert standard array to numpy array:
            roots = np.sort(np.array(roots))

            # Classify roots as ion, unstable, or electron root:
            if len(roots)==1:
                if roots[0]*conversionFactorToEr > 0:
                    root_types=['electron']
                else:
                    root_types=['ion']
            elif len(roots)==3:
                root_types=['ion','unstable','electron']
            else:
                # If the number of roots is not 1 or 3, then flag the case as complicated:
                root_types=['unknown']*len(roots)

                logging.info("Roots found for the radial derivative of the electrostatic potential, for various radial coordinates:")

            logging.info("Root types:"+str(root_types))
         
            # Interpolate the other results onto the ambipolar Er(s):
            logging.info("**** Output quantities at these ambipolar E_r value(s) ****")
            
            output_ambipolar = OrderedDict()
            output_ambipolar['value'] = OrderedDict()
            output_ambipolar['label'] = OrderedDict()
            output_ambipolar['param'] = OrderedDict()
            
            for key in output['value']:
                if whichInterpolator==0:
                    quantityInterpolator = PchipInterpolator(
                        output['value'][ErName]
                        ,output['value'][key])
                else:
                    quantityInterpolator = interp1d(
                        output['value'][ErName]
                        ,output['value'][key])
                value = quantityInterpolator(roots)
                label = output['label'][key]
                
                _output_append_value(output_ambipolar, key, value, label)
                logging.info(key+": "+str(value))
                            
            key = 'root_type'
            label = 'Root type'
            value = root_types
            _output_append_value(output_ambipolar, key, value, label)
            
            
            ambipolaritySolveSucceeded = True
        except Exception as e:
            logging.exception("For some reason, the ambipolarity solve failed. Here is the error message:")
            raise
    else:
        raise Exception("All runs have the same sign of radial current, so we cannot solve for the ambipolar E_r.")

    
    #print(output_ambipolar)
    
    if ambipolaritySolveSucceeded:
        filename = 'ambipolarSolutions.h5'
        filepath = os.path.join(path, filename)
        mirhdf5.dictToHdf5(output_ambipolar, filepath)
        logging.info('Wrote ambipolar solutions to: {}'.format(filepath))
        
        
    return output_ambipolar


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    #logging.basicConfig(level=logging.INFO)
    
    
    run()
    exit(0)
