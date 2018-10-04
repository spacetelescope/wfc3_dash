#!/usr/bin/env python

"""Tests for the ``reduce_dash`` module.

Authors
-------
    Catherine Martlin 2018

Use
---
    These tests can be run via the command line (omit the -s to
    suppress verbose output to stdout):
    ::
        pytest -s test_reduce_dash.py
"""

import pytest

from wfc3tools.wfc3_dash.reduce_dash import *

@pytest.mark.xfail
def test_DashData():
    ''' Testing the init of the class object. These are the current 
    attributes after the init: 
    self.file_name
    self.root
    '''
    testObject = reduce_dash.DashData('insert file path')

    try:
        print(testObject.root)
    except AttributeError:
        print("No `root` attribute.")
    try:
        print(testObject.filename)
    except AttributeError:
        print("No `file_name` attribute.")


def test_split_ima():
    ''' Testing the split_ima function. 
    It creates the following attributes for the class object: 
    self.ima_file
    self.diff
    self.dt
    self.dq
    self.readnoise_2D
    self.file_list
    '''

    try:
        print(testObject.ima_file)
    except AttributeError:
        print("No `ima_file` attribute.")
    try:
        print(testObject.diff)
    except AttributeError:
        print("No `diff` attribute.")
    try:
        print(testObject.dt)
    except AttributeError:
        print("No `dt` attribute.")
    try:
        print(testObject.dq)
    except AttributeError:
        print("No `dq` attribute.")
    try:
        print(testObject.readnoise_2D)
    except AttributeError:
        print("No `readnoise_2D` attribute.")
    try:
        print(testObject.file_list)
    except AttributeError:
        print("No `file_list` attribute.")

    


    assert something == what_you_need_it_to

def test_making_pointing_asn():
    ''' Testing the split_ima function. 
    It creates the following attributes for the class object: 
    self.asn_filename
    '''

    try:
        print(testObject.asn_filename)
    except AttributeError:
        print("No `asn_filename` attribute.")


    assert something == what_you_need_it_to

