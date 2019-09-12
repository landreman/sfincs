"""
A module for dealing with colors.

This module is an extension/modification of parts of:
  matplotlib.colors
  matplotlib.cm
  matplotlib.pyplot

There are a number of things that I don't like about how the matplotlib
versions handle things, so this is my attempt to correct some of these
issues.

Maybe someday I'll try to tackle the matplotlib code directly and create
a replacement for the original code.
"""

import numpy as np
from matplotlib import cm
from matplotlib import colors
from matplotlib.colors import Normalize


defined_gradients = {
    'mirred':{
        'type':'LinearSegmentedColorGradient'
        ,'red':[(0.0, 0.5, 0.5)
                ,(0.5, 1.0, 1.0)
                ,(1.0, 1.0, 1.0)]
        ,'green':[(0.0, 0.0, 0.0)
                  ,(0.5, 0.0, 0.0)
                  ,(1.0, 0.8, 0.8)]
        ,'blue':[(0.0, 0.0, 0.0)
                 ,(0.5, 0.0, 0.0)
                 ,(1.0, 0.8, 0.8)]}
    ,'mirgreen':{
        'type':'LinearSegmentedColorGradient'
        ,'green':[(0.0,1.0,1.0)
                  ,(1.0,0.0,0.0)]}
    ,'mirblue':{
        'type':'LinearSegmentedColorGradient'
        ,'red':[(0.0, 0.0, 0.0)
                ,(0.2, 0.0, 0.0)
                ,(0.4, 0.44, 0.44)
                ,(0.6, 0.44, 0.44)
                ,(0.8, 0, 0.0)
                ,(1.0, 0, 0.0)]
        ,'green':[(0.0, 0.16, 0.16)
                  ,(0.2, 0, 0.0)
                  ,(0.4, 0.6, 0.6)
                  ,(0.6, 0.75, 0.75)
                  ,(0.8, 0.62, 0.62)
                  ,(1.0, 0.41, 0.41)]
        ,'blue':[(0.0, 0.63, 0.63)
                 ,(0.2, 1.0, 1.0)
                 ,(0.4, 0.8, 0.8)
                 ,(0.6, 0.8, 0.8)
                 ,(0.8, 1.0, 1.0)
                 ,(1.0, 0.54, 0.54)]}
    ,'mirRedBlackBlue':{
        'type':'LinearSegmentedColorGradient'
        ,'red':[(0.0, 0.8, 0.8)
                ,(0.5, 0.0, 0.0)
                ,(1.0, 0.0, 0.0)]
        ,'green':[(0.0, 0.0, 0.0)
                  ,(0.5, 0.0, 0.0)
                  ,(1.0, 0.8, 0.8)]
        ,'blue':[(0.0, 0.0, 0.0)
                 ,(0.5, 0.0, 0.0)
                 ,(1.0, 0.8, 0.8)]}
    # mirrainbow is not finished.
    ,'mirrainbow':{
        'type':'LinearSegmentedColorGradient'
        ,'red':[(0.0, 1.0, 1.0)
                ,(0.25, 1.0, 1.0)
                ,(0.5, 0.0, 0.0)
                ,(0.75, 0.0, 0.0)
                ,(1.0, 0.0, 0.0)]
        ,'green':[(0.0, 0.0, 0.0)
                ,(0.25, 0.0, 0.0)
                ,(0.5, 1.0, 1.0)
                ,(0.75, 1.0, 1.0)
                ,(1.0, 0.0, 0.0)]
        ,'blue':[(0.0, 0.0, 0.0)
                ,(0.25, 0.0, 0.0)
                ,(0.5, 0.0, 0.0)
                ,(0.7, 0.0, 0.0)
                ,(1.0, 1.0, 1.0)]}    
}


def getColorGradient(norm=None, cmap=None):
    if not cmap in defined_gradients:

        if not cmap in cm.__dict__.keys():
            raise NameError('Gradient {} not defined'.format(cmap))

        return cm.ScalarMappable(norm, cmap)

    if defined_gradients[cmap]['type'] == 'LinearSegmentedColorGradient':
        gradient = LinearSegmentedColorGradient(norm, defined_gradients[cmap])
    else:
        raise Exception('Gradient {} has wrong type.  Not currently supported'.format(cmap))
        
    return gradient
    
    
class ColorGradient(object):

    def __call__(self, *args):
        return self.to_rgba(*args)

    
class LinearSegmentedColorGradient(ColorGradient):

    rgba_keys = ['red', 'blue', 'green', 'alpha']
            
    def __init__(self, norm=None, segmentdata=None):
        if norm:
            self.norm = norm
        else:
            self.norm = Normalize()
            
        self.setSegmentData(segmentdata)

    
    def _getRgbaComp(self, value, key):

        data = self._segmentdata[key]
        x = data[:, 0]
        y0 = data[:, 1]
        y1 = data[:, 2]

        if value <= x[0]:
            output = y0[0]
        elif value >= x[-1]:
            output = y0[-1]
        else:
            ii = np.searchsorted(x, value)-1

            output = (value - x[ii])/(x[ii+1] - x[ii]) * (y1[ii+1] - y0[ii]) + y0[ii]
            
        return output

    
    def _scalarToRgba(self, value):
        
        return (self._getRgbaComp(self.norm(value), 'red')
                ,self._getRgbaComp(self.norm(value), 'green')
                ,self._getRgbaComp(self.norm(value), 'blue')
                ,self._getRgbaComp(self.norm(value), 'alpha'))

    
    def _arrayToRgba(self, values):
        
        output = np.zeros(values.shape+(4,))

        it = np.nditer(values, flags=['multi_index'])
        while not it.finished:
            output[it.multi_index] = self._scalarToRgba(values[it.multi_index])
            it.iternext()

        return output

    
    def to_rgba(self, value):
               
        if np.isscalar(value):
            output = self._scalarToRgba(value)
        else:
            output = self._arrayToRgba(value)

        return output


    def setSegmentData(self, segmentdata_in):

        self._segmentdata = {}

        for kk in self.rgba_keys:
            if kk in segmentdata_in:
                self._segmentdata[kk] = self._standardizeSegmentData(segmentdata_in[kk])
            elif kk == 'alpha':
                # Set default values for the alpha channel.
                default_data = [[0.0, 1.0, 1.0],[1.0, 1.0, 1.0]]
                self._segmentdata[kk] = self._standardizeSegmentData(default_data)
            else:
                # Set default values for the red, green and blue channels.
                default_data = [[0.0, 0.0, 0.0],[1.0, 0.0, 0.0]]
                self._segmentdata[kk] = self._standardizeSegmentData(default_data)

                
    def _standardizeSegmentData(self, data):

        try:
            array_data = np.array(data)
        except:
            raise TypeError("segment data must be convertable to an array")        

        if not (len(array_data.shape) == 2 and array_data.shape[1] == 3):
            raise ValueError("segment data must be a nx3 array.")

        x = array_data[:, 0]
        if x[0] != 0. or x[-1] != 1.0:
            raise ValueError(
                "segment data points must start with x=0. and end with x=1")
        if np.sometrue(np.sort(x) - x):
            raise ValueError(
                "segment data points must have x in increasing order")

        return array_data
