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
  Main driver for all of the different types of sfincscans.  This is the routine
  that should be called from the command line.
"""

import numpy as np
import sys
import logging

import mirscan_const
import mirscan_util

from mirscan_const import options

def startscan(options=None):
    logging.debug('Entering mirscan.')
    
    mirscan_util.set_defaults_env()
    
    # Load the input file:
    with open(options['input_filename'], 'r') as f:
        input_file = f.readlines()
        
    scan_type = mirscan_util.readScanVariable('scanType', 'int', input_file)
    
    if scan_type == 31:
        import mirscan_31
        mirscan_31.startscan(options=options)
    else:
        raise NotImplementedError('Scantype: {} not yet implemented.'.format(scan_type))
        
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    
    startscan()
    sys.exit(0)

