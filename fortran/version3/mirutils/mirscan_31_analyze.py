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



def start_analysis():
    """
    Loop through each of the Er scans output from mirscan_31 and do an analysis.
    """
   
    # This routines currently assumes the specific mirscan_31 directory structure. 
    
    # Find a list of all completed runs.
    paths_output = glob.glob('**/'+options['output_filename'], recursive=True)
    paths_er = [os.path.dirname(xx) for xx in paths_output]
        
    # If I was really being careful, I would now check that each of these 
    # actually finished by looking inside the hdf5 file.
    
    # Each of these paths should contain an Er scan.
    paths_grid = [os.path.dirname(xx) for xx in paths_er]
    paths_grid = np.unique(paths_grid)
    
    print(paths_grid)
    
    output = OrderedDict()
    output['value'] = OrderedDict()
    output['label'] = OrderedDict()
    output['param'] = OrderedDict()  
    
    for ii, path in enumerate(paths_grid):
        data_er_scan = analyse_er_scan(path)
        
        for key in data_er_scan['value']:
            value = data_er_scan['value'][key]
            label = data_er_scan['label'][key]
            _output_append_value(output, key, value, label)
    
    #for key in output:
    #    print(key)
    #    for key2 in output[key]:
    #        print('  ', key2, output[key][key2])
    
    filename = 'ambipolarSolutions.h5'
    mirhdf5.dictToHdf5(output, filename)
    logging.info('Wrote ambipolar solutions to: {}'.format(filename))    
    
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
            

    
               
def analyse_er_scan(path):
    
    # This routines currently assumes that in the given path are a set of
    # directories that contain the Er scan runs.
    
    filepaths_output = glob.glob(path+'/*/'+options['output_filename'], recursive=True)
    filepaths_input = [os.path.join(os.path.dirname(xx),options['input_filename']) for xx in filepaths_output]
            
    output = OrderedDict()
    output['value'] = OrderedDict()
    output['label'] = OrderedDict()
    output['param'] = OrderedDict()    
        
    for ii in range(len(filepaths_output)):
        
        # This needs to be moved into a separate function.
        # The output from this function is generally useful.
        #
        # Once separated this also should be wrapped in a try/except loop so
        # that we can skip any folders that have errors.
        
        filepath_o = filepaths_output[ii]
        filepath_i = filepaths_input[ii]
        
        obj_file = _open_hdf5(filepath_o)
        try:
            check_run_success(obj_file)
        except:
            logging.exception('Run '+obj_file.filename+' colud not be successfuly read.')
            continue
        

        # This is copied from sfincsScanPlot_2.
        # A lot of work is still needed here.
        with open(filepath_i, 'r') as ff:
            inputfile = ff.readlines()
            
            
        # Extract some information from the input file.
        # Probably this only needs to be done once, but for now it is easier to
        # make this part of the loop.
        inputRadialCoordinate = mirscan_util.readVariable(
            "inputRadialCoordinate"
            ,"int"
            ,inputfile
            ,required=False)
        if inputRadialCoordinate == None:
            inputRadialCoordinate = mirscan_util.readDefault("inputRadialCoordinate", "int", inputfile)
        
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
        
        radius_wish = mirscan_util.readVariable(
            radiusName+"_wish"
            ,"float"
            ,inputfile
            ,required=False)    
        if radius_wish == None:
            radius_wish = mirscan_util.readDefault(radiusName+"_wish","float")

            
        
        # The expression [()] converts from an h5py dataset to a numpy ndarray:
        integerToRepresentTrue = (obj_file["integerToRepresentTrue"][()])
        inputRadialCoordinateForGradients_new = obj_file["inputRadialCoordinateForGradients"][()]
        RHSMode_new = obj_file["RHSMode"][()]
        Nspecies_new = obj_file["Nspecies"][()]
        Zs = obj_file["Zs"][()]
        includePhi1_new = (obj_file["includePhi1"][()] == integerToRepresentTrue)
        if ii == 0:
            inputRadialCoordinateForGradients = inputRadialCoordinateForGradients_new
            RHSMode = RHSMode_new
            Nspecies = Nspecies_new
            includePhi1 = includePhi1_new
        else:
            if inputRadialCoordinateForGradients != inputRadialCoordinateForGradients_new:
                raise Exception("Error! inputRadialCoordinateForGradients is not consistent among runs.")
            if RHSMode != RHSMode_new:
                raise Exception("Error! RHSMode is not consistent among runs.")
            if Nspecies != Nspecies_new:
                raise Exception("Error! Nspecies is not consistent among runs.")
            if includePhi1 != includePhi1_new:
                raise Exception("Error! includePhi1 is not consistent among runs.")
            
        if RHSMode != 1 and RHSMode != 2 and RHSMode != 3:
            raise Exception("Error! sfincsScanPlot is not yet set up for RHSMode = "+str(RHSMode))          

        doAmbipolaritySolve = (RHSMode == 1)          

        key = radiusName+'_wish'
        label = 'Radius Wish'
        value = radius_wish
        _output_append_value(output, key, value, label)
        
        key = radiusName
        label = 'Radius'
        value = obj_file[key][()]
        _output_append_value(output, key, value, label)
  
    
        if inputRadialCoordinateForGradients==0:
            varName = "dPhiHatdpsiHat"
        elif inputRadialCoordinateForGradients==1:
            varName = "dPhiHatdpsiN"
        elif inputRadialCoordinateForGradients==2:
            varName = "dPhiHatdrHat"
        elif inputRadialCoordinateForGradients==3:
            varName = "dPhiHatdrN"
        elif inputRadialCoordinateForGradients==4:
            varName = "Er"
        else:
            raise Exception("Error! Invalid inputRadialCoordinateForGradients.")
        
        key = varName
        label = varName
        value = obj_file[varName][()]
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
    
    #if not atLeastOneDirectorySucceeded:
    #   raise("Error! There do not seem to be any completed sfincs jobs in subdirectories of this directory.")
    
    conversionFactorToddpsiHat = obj_file['dPhiHatdpsiHat'][()] / obj_file[varName][()]
    conversionFactorToddpsiN   = obj_file['dPhiHatdpsiN'][()]   / obj_file[varName][()]
    conversionFactorToddrHat   = obj_file['dPhiHatdrHat'][()]   / obj_file[varName][()]
    conversionFactorToddrN     = obj_file['dPhiHatdrN'][()]     / obj_file[varName][()]
    conversionFactorToEr       = obj_file['Er'][()]             / obj_file[varName][()]
    
    
    # Sort by Er:
    sort_index = np.argsort(output['value'][varName])
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
    if len(output['value'][varName]) < 3:
        whichInterpolator = 1
    else:
        whichInterpolator = 0
    
    ambipolaritySolveSucceeded = False
    if doAmbipolaritySolve:
        if max(output['value']['radialCurrent']) > 0 and min(output['value']['radialCurrent']) < 0:
            try:
                # We have runs that straddle j_psi = 0, so we can find at least 1 ambipolar solution.
                logging.info("Solving for E_r.")
                
                # Note: For interpolation it does not work well to use scipy's "interp1d" function, 
                #       since this tends to give large overshoots for the typical data we get in an 
                #       Er scan (with sharp structure near Er=0).
                if whichInterpolator==0:
                    interpolator = PchipInterpolator(output['value'][varName], output['value']['radialCurrent'])
                else:
                    interpolator = interp1d(output['value'][varName], output['value']['radialCurrent'])
                 
                # In the next line, 500 is arbitrary, and could be replaced by another large integer.
                Er_fine = np.linspace(min(output['value'][varName]), max(output['value'][varName]), num=500)
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
                message = 'd PhiHat / d psiHat = '+' '.join(
                    ['{}'.format(xx) for xx in (roots * conversionFactorToddpsiHat)])
                logging.info(message)
                message = 'd PhiHat / d psiN   = '+' '.join(
                    ['{}'.format(xx) for xx in (roots * conversionFactorToddpsiN)])
                logging.info(message)
                message = 'd PhiHat / d rHat   = '+' '.join(
                    ['{}'.format(xx) for xx in (roots * conversionFactorToddrHat)])
                logging.info(message)
                message = 'd PhiHat / d rN     = '+' '.join(
                    ['{}'.format(xx) for xx in (roots * conversionFactorToddrN)])
                logging.info(message)
                message = 'Er                  = '+' '.join(
                    ['{}'.format(xx) for xx in (roots * conversionFactorToEr)])
                logging.info(message)
    
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
                            output['value'][varName]
                            ,output['value'][key])
                    else:
                        quantityInterpolator = interp1d(
                            output['value'][varName]
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
    else:
        logging.info("Not doing an ambipolarity solve since RHSMode > 1.")

    
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
    
    
    start_analysis()
    exit(0)
