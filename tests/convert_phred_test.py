#!/usr/bin/env python

import sys
sys.path.append("..")

import bioinfo

def test_convert_phred():
    assert bioinfo.convert_phred("I") == 40, "wrong phred score for 'I'"
    assert bioinfo.convert_phred("C") == 34, "wrong phred score for 'C'"
    assert bioinfo.convert_phred("2") == 17, "wrong phred score for '2'"
    assert bioinfo.convert_phred("@") == 31, "wrong phred score for '@'"
    assert bioinfo.convert_phred("$") == 3, "wrong phred score for '$'"
