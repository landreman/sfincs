#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
"""
A set of tools to deal with HDF5 files.

Author:
  Novimir Antoniuk Pablant
    - npablant@pppl.gov
    - novimir.pablant@amicitas.com

Important notes:
  Python specific types/objects cannot in general be directly translated to 
  HDF5 data types.  Generally the use of non compatable objects/types will
  produce an error.  There are however a few exceptions that have been added
  to mirHDF5 for convenience.

  None
      The python None object will be saved as False, with specific 
      HDF5 attributes to allow proper restortation in python.

  Dict
      Python dictionaries are saved using HDF5 group and data entries.

  List
      Python lists are saved using HDF5 group and data entries.  Entry
      names are automatically generated, and an HDF5 attribute is set
      to allow the python list to be restored.

ToDo:
  Currently unicode key names are not fully supported, at least for
  python-2.7.  I currently get a TypeError warning when trying to
  save the key order.  I need to look into this a bit to see if simply
  converting to str will fix this.
"""

# For python 2.7/3.0 compatability.
from __future__ import print_function
from __future__ import unicode_literals
str_builtin = str
str = ''.__class__

import logging
import pathlib
import zipfile
import re

import numpy as np
import h5py
from collections import OrderedDict

def dictToHdf5(dick, filename=None, driver=None, **kwargs):
    """
    Save a python dictionary as an HDF5 data file.

    This will only work for dictionaries where the keys are all strings. The
    order of the dictionary will always be preserved.

    programming notes
    -----------------

    In order to retain the order of the input dictionary I create an dictionary
    order attribute as part of the HDF5 root '/'. The h5py library will by default 
    always return the groups before the datasets, so I need to save the original
    order somehow.
    """
    
    # I need to check if the file exists here.  
    # This can also be done by using the correct flags in h5py.File.
    file = h5py.File(filename, 'w', driver=driver)

    if not isinstance(dick, dict):
        raise Exception('Incorrect input type. Dictionary expected.')

    _addDictToHdf5(file, dick, **kwargs)

    file.close()


def hdf5ToDict(filename, driver=None, include=None, exclude=None):
    """
    Read an HDF5 file into a dictionary.

    This will only work for simple HDF files in which we only care
    about the data name and value.

    keywords
    --------

    exclude
      String or list of strings. Strings must be regular expressions.
      Anykeys matching an exclude pattern will be excluded in the
      returned dictionary. 

    include
      String or list of strings. Strings must be regular expressions.
      Any excluded keys matching an include pattern will be included
      in the returned dictionary. If no exclude pattern is given then
      exclude will be set to '.*'. The order of evaluation is exclude
      patterns first, then include patterns.

    WARNING:
      This routine has some unicode issues in python 2.  In python 2
      the dictionary keys  will always be returned in unicode regardless
      of the original keys.
    """

    include, exclude = _setupIncludeExclude(include, exclude)
    
    file = h5py.File(filename, 'r', driver=driver)

    dick = _createNewItemFromHdf5(file, include=include, exclude=exclude)
    
    file.close()

    return dick


def dictToHdf5Zip(dick, filename=None, **kwargs):
    """
    Write a dictionary into a zipped hdf5 file.
    
    This works by:
      1. writing a hdf5 file.
      2. moving that file into a zip archive.
      3. delete the hdf5 file.

    This will overwrite and delete existing files without any warning!
    """
    
    p = pathlib.Path(filename)
    if p.suffix == '.zip':
        file_z = filename
        file_h = str(p.parent/p.stem)
    else:
        file_z = filename+'.zip'
        file_h = filename

    # is_file_h = pathlib.Path(path_h).exists()
    # is_file_z = pathlib.Path(path_z).exists()
    
    dictToHdf5(dick, filename=file_h, **kwargs)

    with zipfile.ZipFile(
            file_z
            ,mode='w'
            ,compression=zipfile.ZIP_DEFLATED
            ,compresslevel=6
            ) as ff_z:
        ff_z.write(file_h, pathlib.Path(file_h).name)

    pathlib.Path(file_h).unlink()

    
def hdf5ZipToDict(filename, **kwargs):
    """
    This works by:
      1. Extract hdf5 file from zip archive.
      2. Read hdf5 file.
      3. delete the hdf5 file.

    This will overwrite and delete existing files without any warning!
    """

    p = pathlib.Path(filename)
    if p.suffix == '.zip':
        file_z = filename
        file_h = str(p.parent/p.stem)
    else:
        file_z = filename+'.zip'
        file_h = filename

    # is_file_h = pathlib.Path(path_h).exists()
    # is_file_z = pathlib.Path(path_h).exists()

    # Extract the hdf5 file.
    # Here we assume that there is only one file in the zip archive.
    with zipfile.ZipFile(file_z, mode='r') as ff_z:
        # Strip out any path information.
        info = ff_z.infolist()[0]
        info.filename = pathlib.Path(info.filename).name
        ff_z.extract(ff_z.infolist()[0], p.parent)
        
    dick = hdf5ToDict(file_h, **kwargs)
    pathlib.Path(file_h).unlink()

    return dick

    
def _addItemToHdf5(group, key, item, compression=None, compression_opts=None):
    if compression is True: compression = 'lzf'

    try:
        if isinstance(item, dict):
            new_group = group.create_group(key)
            _addDictToHdf5(
                new_group
                ,item
                ,compression=compression
                ,compression_opts=compression_opts)
        elif isinstance(item, list):
            new_group = group.create_group(key)
            _addListToHdf5(
                new_group
                ,item
                ,compression=compression
                ,compression_opts=compression_opts)
        elif item is None:
            group[key] = False
            group[key].attrs['_mirhdf5 python None'] = True
        else:
            if np.isscalar(item):
                group.create_dataset(
                    key
                    ,data=item)
            else:
                group.create_dataset(
                    key
                    ,data=item
                    ,compression=compression
                    ,compression_opts=compression_opts)
    except TypeError:
        logging.exception('Could not add key "{}" of type {} to hdf5 file.'.format(key, type(item)))


def _addDictToHdf5(group, dick, **kwargs):

    keys = [key.encode() for key in dick.keys()]
    group.attrs['_mirhdf5 python object type'] = 'OrderedDict'.encode()
    try:
        group.attrs['_mirhdf5 dictionary order'] = keys
    except TypeError:
        logging.error('Could not save dictionary key order. keys of unsupported data type.')


    for key in dick:
        _addItemToHdf5(group, key, dick[key], **kwargs)


def _addListToHdf5(group, in_list, **kwargs):

    group.attrs['_mirhdf5 python object type'] = 'list'.encode()

    for ii, item in enumerate(in_list):
        key = '{:04d}'.format(ii)
        _addItemToHdf5(group, key, item, **kwargs)


def _addHdf5ToDict(
        group
        ,dick
        ,attrs=None
        ,group_type=None
        ,include=None
        ,exclude=None
        ):

    if group_type != 'OrderedDict':
        raise Exception('HDF5 parsing error. Expected a OrderedDict type, got {} instead.'.format(group.attrs['_mirhdf5 python object type']))

    if '_mirhdf5 dictionary order' in attrs:
        key_list = attrs['_mirhdf5 dictionary order']
    else:
        key_list = group.keys()

    for key in key_list:
        # Convert bytes key to approprate Python 2/3 string type.
        key_string = key
        if hasattr(key, 'decode'):
            key_string = key.decode()
                
        if _keyIncluded(key_string, include, exclude):
            new_item = _createNewItemFromHdf5(group, key, include, exclude)
            dick[key_string] = new_item


def _addHdf5ToList(
        group
        ,in_list
        ,attrs=None
        ,group_type=None
        ,include=None
        ,exclude=None
        ):

    if group_type != 'list':
        raise Exception('HDF5 parsing error. Expected a list type.')
    
    keys = group.keys()
    for key in keys:
        new_item = _createNewItemFromHdf5(group, key, include, exclude)
        in_list.append(new_item)


def _createNewItemFromHdf5(group, key=None, include=None, exclude=None):

    if key is None:
        data = group
    else:
        data = group[key]
    attrs = data.attrs
    
    if isinstance(data, h5py.Group):
        if '_mirhdf5 python object type' in attrs:
            group_type = attrs['_mirhdf5 python object type']
            if hasattr(group_type, 'decode'):
                group_type = group_type.decode()
        else:
            group_type = 'OrderedDict'

        if group_type == 'OrderedDict':
            new_item = OrderedDict()
            _addHdf5ToDict(data, new_item, attrs, group_type, include, exclude)
        elif group_type == 'list':
            new_item = list()
            _addHdf5ToList(data, new_item, attrs, group_type, include, exclude)
        else:
            raise Exception('Unknown group type: {}'.format(group_type))
    elif data is False:
        new_item = False
        if '_mirhdf5 python None' in attrs:
            if attrs['_mirhdf5 python None'] == True:
                new_item = None
    else:
        new_item = data[()]

    return new_item


def _setupIncludeExclude(include, exclude):
        
    # If a single string was given, turn it into a list.
    # Then compile all the regular expressions.
    if include is not None:
        if isinstance(include, str):
            include = [include]
        include = [re.compile(xx) for xx in include]

        # if only include was provide, then exclude everything else.
        if exclude is None:
            exclude = ['.*']
            
    if exclude is not None:
        if isinstance(exclude, str):
            exclude = [exclude]
        exclude = [re.compile(xx) for xx in exclude]
        
    return include, exclude


def _keyIncluded(key, include, exclude):

    status = True
    if exclude:
        for regex in exclude:
            if regex.match(key):
                status = False
    if include:
        for regex in include:
            if regex.match(key):
                status = True
                
    return status
