#!/usr/bin/env python 

DNA_bases = "ATCGatcg"
RNA_bases = "AUCGaucg"

DNA_comp = {'A':'T', 'a':'t',
            'T':'A', 't':'a',
            'C':'G', 'c':'g',
            'G':'C', 'g':'c',
            'N':'N', 'n':'n'}
RNA_comp = {'A':'U', 'a':'u',
            'U':'A', 'u':'a',
            'C':'G', 'c':'g',
            'G':'C', 'g':'c',
            'N':'N', 'n':'n'}

codons = {
        'AAA': ['Lys', 'K',	'Lysine'],
        'AAC': ['Asn', 'N', 'Asparagine'],
        'AAG': ['Lys', 'K', 'Lysine'],
        'AAT': ['Asn', 'N', 'Asparagine'],
        'ACA': ['Thr', 'T', 'Threonine'],
        'ACC': ['Thr', 'T', 'Threonine'],
        'ACG': ['Thr', 'T', 'Threonine'],
        'ACT': ['Thr', 'T', 'Threonine'],
        'AGA': ['Arg', 'R', 'Arginine'],
        'AGC': ['Ser', 'S', 'Serine'],
        'AGG': ['Arg', 'R', 'Arginine'],
        'AGT': ['Ser', 'S', 'Serine'],
        'ATA': ['Ile', 'I', 'Isoleucine'],
        'ATC': ['Ile', 'I', 'Isoleucine'],
        'ATG': ['Met', 'M', 'Methionine'],
        'ATT': ['Ile', 'I', 'Isoleucine'],
        'CAA': ['Gln', 'Q', 'Glutamine'],
        'CAC': ['His', 'H', 'Histidine'],
        'CAG': ['Gln', 'Q', 'Glutamine'],
        'CAT': ['His', 'H', 'Histidine'],
        'CCA': ['Pro', 'P', 'Proline'],
        'CCC': ['Pro', 'P', 'Proline'],
        'CCG': ['Pro', 'P', 'Proline'],
        'CCT': ['Pro', 'P', 'Proline'],
        'CGA': ['Arg', 'R', 'Arginine'],
        'CGC': ['Arg', 'R', 'Arginine'],
        'CGG': ['Arg', 'R', 'Arginine'],
        'CGT': ['Arg', 'R', 'Arginine'],
        'CTA': ['Leu', 'L', 'Leucine'],
        'CTC': ['Leu', 'L', 'Leucine'],
        'CTG': ['Leu', 'L', 'Leucine'],
        'CTT': ['Leu', 'L', 'Leucine'],
        'GAA': ['Glu', 'E', 'Glutamic_acid'],
        'GAC': ['Asp', 'D', 'Aspartic_acid'],
        'GAG': ['Glu', 'E', 'Glutamic_acid'],
        'GAT': ['Asp', 'D', 'Aspartic_acid'],
        'GCA': ['Ala', 'A', 'Alanine'],
        'GCC': ['Ala', 'A', 'Alanine'],
        'GCG': ['Ala', 'A', 'Alanine'],
        'GCT': ['Ala', 'A', 'Alanine'],
        'GGA': ['Gly', 'G', 'Glycine'],
        'GGC': ['Gly', 'G', 'Glycine'],
        'GGG': ['Gly', 'G', 'Glycine'],
        'GGT': ['Gly', 'G', 'Glycine'],
        'GTA': ['Val', 'V', 'Valine'],
        'GTC': ['Val', 'V', 'Valine'],
        'GTG': ['Val', 'V', 'Valine'],
        'GTT': ['Val', 'V', 'Valine'],
        'TAA': ['Stp', 'O', 'Stop'],
        'TAC': ['Tyr', 'Y', 'Tyrosine'],
        'TAG': ['Stp', 'O', 'Stop'],
        'TAT': ['Tyr', 'Y', 'Tyrosine'],
        'TCA': ['Ser', 'S', 'Serine'],
        'TCC': ['Ser', 'S', 'Serine'],
        'TCG': ['Ser', 'S', 'Serine'],
        'TCT': ['Ser', 'S', 'Serine'],
        'TGA': ['Stp', 'O', 'Stop'],
        'TGC': ['Cys', 'C', 'Cysteine'],
        'TGG': ['Trp', 'W', 'Tryptophan'],
        'TGT': ['Cys', 'C', 'Cysteine'],
        'TTA': ['Leu', 'L', 'Leucine'],
        'TTC': ['Phe', 'F', 'Phenylalanine'],
        'TTG': ['Leu', 'L', 'Leucine'],
        'TTT': ['Phe', 'F', 'Phenylalanine']}

class sequence(): 
    """
    Nucleic acid sequence, default DNA, may be RNA.
    If RNA, user must define. 
    No 'N' bases allowed.
    """ 

    DEFAULT_KIND = "DNA"

    def __init__(self, seq, kind = DEFAULT_KIND):
        self.kind = kind
        self.seq = seq
        assert self.validateseq()
        self.rc = self.rev_comp()

    def validateseq(self): 
        if self.kind == "DNA": return set(DNA_bases).issuperset(self.seq)
        if self.kind == "RNA": return set(RNA_bases).issuperset(self.seq)

    def gc_content(self):
        length: int = len(self.seq)
        gcount: int = self.seq.count('G') + self.seq.count('g')
        ccount: int = self.seq.count('C') + self.seq.count('c')
        acount: int = self.seq.count('A') + self.seq.count('a')
        tcount: int = self.seq.count('T') + self.seq.count('t')
        return (ccount + gcount)/length

    def rev_comp(self):
        '''Returns the reverse complement of seq.'''

        reversecomp: str = ""
        if self.kind == "DNA":
            for index in range(len(self.seq)):
                reversecomp = reversecomp + DNA_comp[self.seq[-1*(index + 1)]]
        if self.kind == "RNA": 
            for index in range(len(self.seq)):
                reversecomp = reversecomp + RNA_comp[self.seq[-1*(index + 1)]]
        return reversecomp

    def whatami(self):
        print(self.kind) 

    def uppercase(self): 
        self.seq = self.seq.upper()
        pass

if __name__ == "__main__":
    tseq = sequence("ATGGGCGC")
    tseq.whatami()
    tseq2 = sequence("AUGGGCGC", "RNA") 
    print(tseq2.rc)
    tseq3 = sequence("ATCGatcgg") 
    tseq3.uppercase()
    print(tseq3.seq)
    print(tseq3.gc_content())
