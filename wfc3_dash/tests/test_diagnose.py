#!/usr/bin/env python

"""Tests for the ``diagnose`` module.

Authors
-------
    Catherine Martlin 2018
    
Use
---
    These tests can be run via the command line (omit the -s to
    suppress verbose output to stdout):
    ::
        pytest -s test_diasnose.py
"""

import pytest

from wfc3tools.wfc3_dash.diagnose import function

@pytest.mark.xfail
def test_function():
    ''' Do something that is asserting to test things.'''


    assert something == what_you_need_it_to

