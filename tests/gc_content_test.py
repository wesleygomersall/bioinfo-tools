#!/usr/bin/env python

import sys
sys.path.append("..")

import bioinfo

def test_gc_content():
    assert bioinfo.gc_content("GCGCGC") == 1
    assert bioinfo.gc_content("AATTATA") == 0
    assert bioinfo.gc_content("GCATCGAT") == 0.5
