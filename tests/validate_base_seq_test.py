#!/usr/bin/env python

import sys
sys.path.append("..")

import bioinfo

def test_validate_base_seq():
    assert bioinfo.validate_base_seq("AATAGAT"), "Validate base seq does not work on DNA"
    assert bioinfo.validate_base_seq("AAUAGAU", True), "Validate base seq does not work on RNA"
    assert bioinfo.validate_base_seq("R is the best!")==False, "Not a DNA string"
    assert bioinfo.validate_base_seq("aatagat"), "Validate base seq does not work on lowercase DNA"
    assert bioinfo.validate_base_seq("aauagau", True), "Validate base seq does not work on lowercase RNA"
    assert bioinfo.validate_base_seq("TTTTtttttTTT")
