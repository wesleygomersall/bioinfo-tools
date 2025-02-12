#!/usr/bin/env python

import sys
sys.path.append("..")

import bioinfo

def test_calc_median():
    assert bioinfo.calc_median([1,2,3]) == 2
    assert bioinfo.calc_median([5,6,7,8]) == 6.5
    assert bioinfo.calc_median([1,1,1,1,1,1,1,1,100]) == 1
    assert bioinfo.calc_median([7]) == 7
    assert bioinfo.calc_median([50,100]) == 75
