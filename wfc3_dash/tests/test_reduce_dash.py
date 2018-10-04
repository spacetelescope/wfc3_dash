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

from wfc3tools.wfc3_dash.reduce_dash import function

@pytest.mark.xfail
def test_function():
    ''' Do something that is asserting to test things.'''


    assert something == what_you_need_it_to

