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
import logging
import os
from collections import OrderedDict

import mirhdf5
import mirplot


def options_default():
    options = OrderedDict()

    options['root'] = None
    
    options['xbound'] = [-0.1, 1.1]
    
    options['ybound_n'] = None
    options['ybound_T'] = None
    options['ybound_Er'] = None
    
    return options


def create_root_mask(data, root=None, keepone=False):
    rN = np.array(data['value']['rN']).round(6)

    if root is not None:
        mask = np.zeros(len(rN), dtype=bool)
        rho_unique = np.unique(rN)
        for ii, rho_dummy in enumerate(rho_unique):
            temp = np.nonzero(rN == rho_dummy)[0]
            if (len(temp) > 1) or (not keepone):
                w = (data['value']['root_type'][temp] == root)
                mask[temp[w]] = True
            else:
                mask[temp] = True
    else:
        mask = np.ones(len(rN), dtype=bool)
    
    return mask


def start_plot(path=None, options=None):
    """
    Plot the results from a mirscan_31.
    """

    if path is None: path = ''

    opt = options_default()
    if options is not None:
        opt.update(options)        
    
    filename = 'ambipolarSolutions.h5'
    filepath = os.path.join(path, filename)
    output = mirhdf5.hdf5ToDict(filepath)
    logging.info('Reading ambipolar solutions from: {}'.format('ambipolarSolutions.h5'))
    
    
    color = []
    for root in output['value']['root_type']:
        if root == 'ion':
            color.append('blue')
        elif root == 'electron':
            color.append('red')
        else:
            color.append('black')


    mask = create_root_mask(output, opt['root'])
    
    obj_plot = mirplot.MirPlot()
    obj_plot.append({
        'name':'n'
        ,'type':'line'
        ,'x':output['value']['rN']
        ,'y':output['value']['nHats_1']
        ,'xlabel':'rho'
        ,'ylabel':'n, T'
        ,'label':'nHats_1'
        ,'legend':True
        ,'xbound':opt['xbound']
        ,'ybound':opt['ybound_n']
        ,'color':'black'
        })
    obj_plot.append({
        'name':'n'
        ,'type':'hline'
        ,'y':[0]
        ,'color':'black'
        ,'alpha':0.2
        ,'linestyle':'--'
        })
    
    obj_plot.append({
        'name':'T'
        ,'type':'line'
        ,'x':output['value']['rN']
        ,'y':output['value']['THats_1']
        ,'xlabel':'rho'
        ,'label':'THats_1'
        ,'xbound':opt['xbound']
        ,'ybound':opt['ybound_T']
        ,'color':'darkred'
        })
    obj_plot.append({
        'name':'T'
        ,'type':'line'
        ,'x':output['value']['rN']
        ,'y':output['value']['THats_2']
        ,'xlabel':'rho'
        ,'label':'THats_2'
        ,'color':'darkblue'
        })
    obj_plot.append({
        'name':'T'
        ,'type':'hline'
        ,'y':[0]
        ,'color':'black'
        ,'alpha':0.2
        ,'linestyle':'--'
        })
    
    obj_plot.append({
        'name':'Er'
        ,'type':'scatter'
        ,'x':output['value']['rN']
        ,'y':output['value']['Er']
        ,'color':color
        ,'markersize':50.0
        ,'xlabel':'rho'
        ,'ylabel':'Er'
        ,'xbound':opt['xbound']
        ,'ybound':opt['ybound_Er']
        })
    obj_plot.append({
        'name':'Er'
        ,'type':'hline'
        ,'y':[0]
        ,'color':'black'
        ,'alpha':0.2
        ,'linestyle':'--'
        })
    obj_plot.setFigureProperties({'figure_title':path})
    obj_plot.plotToScreen()

    return obj_plot
    
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    #logging.basicConfig(level=logging.INFO)
    
    
    start_plot()
    exit(0)
