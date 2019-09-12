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

import mirhdf5
import mirplot

def start_plot():
    """
    Plot the results from a mirscan_31.
    """
    
    filename = 'ambipolarSolutions.h5'
    output = mirhdf5.hdf5ToDict(filename)
    logging.info('Reading ambipolar solutions from: {}'.format('ambipolarSolutions.h5'))
    
    
    color = []
    for root in output['value']['root_type']:
        if root == 'ion':
            color.append('blue')
        elif root == 'electron':
            color.append('red')
        else:
            color.append('black')
    
    xbound = [-0.1, 1.1]
    
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
        ,'xbound':xbound
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
        ,'xbound':xbound
        ,'color':'darkred'
        })
    obj_plot.append({
        'name':'T'
        ,'type':'line'
        ,'x':output['value']['rN']
        ,'y':output['value']['THats_2']
        ,'xlabel':'rho'
        ,'label':'THats_2'
        ,'xbound':xbound
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
        ,'xbound':xbound
        })
    obj_plot.append({
        'name':'Er'
        ,'type':'hline'
        ,'y':[0]
        ,'color':'black'
        ,'alpha':0.2
        ,'linestyle':'--'
        })
    obj_plot.plotToScreen()
    
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    #logging.basicConfig(level=logging.INFO)
    
    
    start_plot()
    exit(0)
