#!/usr/bin/env python

import sys
sys.path.append("..")

import bioinfo

def test_qual_score():
    assert bioinfo.qual_score("A") == 32.0, "wrong average phred score for 'A'"
    assert bioinfo.qual_score("AC") == 33.0, "wrong average phred score for 'AC'"
    assert bioinfo.qual_score("@@##") == 16.5, "wrong average phred score for '@@##'"
    assert bioinfo.qual_score("EEEEAAA!") == 30.0, "wrong average phred score for 'EEEEAAA!'"
    assert bioinfo.qual_score("$") == 3.0, "wrong average phred score for '$'"
