#!/usr/bin/env python

import os
import argparse 
import bioinfo


def get_args():
    parser = argparse.ArgumentParser(description=" Program for translation of mRNA sequences. Outputs fasta format amino acid translation (to standard output) of input string or fasta file of sequences. ")
    parser.add_argument("-s", "--sequence", help=" mRNA sequence to translate to amino acid sequence. This can be either a fasta file or a sequence string. ", type=str, required=True) 
    parser.add_argument("-e", "--end", help=" Use if translation should ignore stop codons in sequence. Will be represented by an 'O' in amino acid sequence. ", action='store_true')
    return parser.parse_args()

codon_table =  {'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAT':'N',
                'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AGA':'R', 'AGC':'S', 'AGG':'R', 'AGT':'S',
                'ATA':'I', 'ATC':'I', 'ATG':'M', 'ATT':'I',
                'CAA':'Q', 'CAC':'H', 'CAG':'Q', 'CAT':'H',
                'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                'GAA':'E', 'GAC':'D', 'GAG':'E', 'GAT':'D',
                'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                'TAA':'O', 'TAC':'Y', 'TAG':'O', 'TAT':'Y', # stop codon gives letter 'O'
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                'TGA':'O', 'TGC':'C', 'TGG':'W', 'TGT':'C',
                'TTA':'L', 'TTC':'F', 'TTG':'L', 'TTT':'F'}

def reading_frame(mRNA_seq: str) -> list: 
    """ Given an mRNA sequence, finds reading frame for the sequence by looking for the start codon 'ATG', returning a list of the 3 base codons. """
    frame = list() 
    reading = False
    base1, base2, base3 = 'N', 'N', 'N'
    index = 0
    inframe = 0
    while True: 
        
        base1 = base2
        base2 = base3
        base3 = mRNA_seq[index]
        current_codon = base1 + base2 + base3

        if reading:
            inframe += 1
            if inframe == 3: 
                frame.append(current_codon) 
                inframe = 0
        
        if not reading and current_codon.upper() in ["ATG", "AUG"]:
            reading = True
            frame.append(current_codon) 
        
        index += 1
        if index == len(mRNA_seq): 
            break

    return(frame) 

def translate_codons(codons: list, stop: bool = True) -> str: 
    """ Given a list of nucleic acid codons, will return a string of the amino acid sequence produced by these codons. 
    Codons can have either T or U in sequence, this function treats these bases the same. 
    Pass an optional boolean for whether or not to stop at the first stop codon. Default is true.
    For 'False', the stop codon will be "translated" to 'O'. """
    aasequence: str = ""
    for codon in codons: 
        aa = codon_table[codon.upper().replace("U", "T")]
        if stop and aa == 'O':
            break
        aasequence = f"{aasequence}{aa}"
    return aasequence

if __name__ == "__main__": 
    args = get_args()

    if os.path.exists(args.sequence): 
        with open(args.sequence, 'r') as fin:
            read: int = ""
            while True:
                line: str = fin.readline()
                if line == "":
                    translation = translate_codons(reading_frame(read))
                    print(translation)
                    break
                if line[0] == '>':
                    if read != "":
                        translation = translate_codons(reading_frame(read))
                        print(translation)
                    print(f"{line.strip('\n')}_aatranslated")
                    read: str = ""
                else:
                    read = read + line.strip().strip()
    
    else:
        seq = args.sequence
        if bioinfo.validate_base_seq(seq, False, True) or bioinfo.validate_base_seq(seq, True, True):
            protseq = translate_codons(reading_frame(seq))
            print(protseq)
        else: 
            print("{args.sequence} is not a valid DNA or RNA sequence")
