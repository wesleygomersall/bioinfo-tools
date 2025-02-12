#!/usr/bin/env python

import sys
sys.path.append("..")

from translate_mRNA import reading_frame
from translate_mRNA import translate_codons

if __name__ == "__main__": 
    assert reading_frame("NNNNATGAAATTTCCCGGG") == ["ATG", "AAA", "TTT", "CCC", "GGG"]
    assert reading_frame("ATGAAATTTCCCGGG") == ["ATG", "AAA", "TTT", "CCC", "GGG"]
    assert reading_frame("NNNNNCTGTGCTGAGTNNNNATGAAATTTCCCGGG") == ["ATG", "AAA", "TTT", "CCC", "GGG"]
    assert reading_frame("NNATGCTGACGGGATTG") == ["ATG", "CTG", "ACG", "GGA", "TTG"]

    assert translate_codons(["ATG", "AAA", "TTT", "CCC", "GGG"]) == "MKFPG"
