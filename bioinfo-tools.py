#!/usr/bin/env python 

DNA_BASES = "ATCGatcg"

RNA_BASES = "AUCGaucg"

DNA_COMP = {'A':'T', 'a':'t',
            'T':'A', 't':'a',
            'C':'G', 'c':'g',
            'G':'C', 'g':'c',
            'N':'N', 'n':'n'}

RNA_COMP = {'A':'U', 'a':'u',
            'U':'A', 'u':'a',
            'C':'G', 'c':'g',
            'G':'C', 'g':'c',
            'N':'N', 'n':'n'}

CODONS = {'AAA': ['Lys', 'K',	'Lysine'],
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

IUPAC_REGEX = {'A': '[Aa]', 'a': '[Aa]',
               'C': '[Cc]', 'c': '[Cc]',
               'T': '[Tt]', 't': '[Tt]', 
               'G': '[Gg]', 'g': '[Gg]',
               'U': '[UuTt]', 'u': '[UuTt]',
               'W': '[AaTt]', 'w': '[AaTt]',
               'S': '[CcGg]', 's': '[CcGg]',
               'M': '[AaCc]', 'm': '[AaCc]',
               'K': '[GgTt]', 'k': '[GgTt]',
               'R': '[AaGg]', 'r': '[AaGg]',
               'Y': '[CcTt]', 'y': '[CcTt]',
               'B': '[CcGgTt]', 'b': '[CcGgTt]',
               'D': '[AaGgTt]', 'd': '[AaGgTt]',
               'H': '[AaCcTt]', 'h': '[AaCcTt]',
               'V': '[AaCcGg]', 'v': '[AaCcGg]',
               'N': '[AaCcTtUuGg]', 'n': '[AaCcTtUuGg]'}

class NucleicAcid(): 
    '''Nucleic acid sequence, default DNA.
    If RNA, user must define is_DNA as False. 
    No undefined bases allowed.

    Attributes:
    is_dna: (bool) Sequence is DNA (RNA else)
    seq:    (str) Sequence
    length: (int) Sequence length
    ''' 

    def __init__(self, sequence: str, is_DNA: bool = True):
        self.is_dna: str = is_DNA
        self.seq: str = sequence.upper()
        self.length: int = len(sequence)
        assert self.validateseq()

    def validateseq(self): 
        '''Must be a sequence containing only A T C or G if DNA,
        or a A U C or G if RNA.
        Case insensitive.'''
        if self.is_dna == True: return set(DNA_BASES).issuperset(self.seq)
        if self.is_dna == False: return set(RNA_BASES).issuperset(self.seq)

    def gc_content(self) -> float:
        '''Calculate GC content of sequence.'''
        gcount: int = self.seq.count('G') + self.seq.count('g')
        ccount: int = self.seq.count('C') + self.seq.count('c')
        return (ccount + gcount) / self.length

    def rev_comp(self) -> str:
        '''Returns the reverse complement of seq.'''
        reversecomp: str = ""
        if self.is_dna == True:
            for index in range(len(self.seq)):
                reversecomp = reversecomp + DNA_COMP[self.seq[-1*(index + 1)]]
        if self.is_dna == False: 
            for index in range(len(self.seq)):
                reversecomp = reversecomp + RNA_COMP[self.seq[-1*(index + 1)]]
        return reversecomp

    def transcribe(self):
        '''Transcribes DNA sequence into RNA in-place. 
        Assumes DNA sequence is coding.
        Changes is_dna to False'''

        assert self.is_dna == True 
        rna = ""
        for i in range(len(self.seq)): 
            if self.seq[i] == "T": rna = rna + "U"
            else: rna = rna + self.seq[i]
        self.seq = rna
        self.is_dna = False 

    def rev_transcribe(self):
        '''Reverse-transcribes RNA sequence into DNA in-place. 
        Changes is_DNA to True'''

        assert self.is_dna == False 
        dna = ""
        for i in range(len(self.seq)): 
            if self.seq[i] == "U": dna = dna + "T"
            else: dna = dna + self.seq[i]
        self.seq = dna
        self.is_dna = True 

    def translate(self) -> str: 
        '''Translate DNA nuc_acid. 
        Sequence must have a start codon to find the reading frame.
        Translates up to first stop codon encountered in seq.'''

        assert self.is_dna ==True 
        
        # to begin: empty protein str and not translating
        protein: str = ""; translating: bool = False 
        frame = self.seq[0:3]; i = 3

        if frame == 'ATG': 
            translating = True; protein = "M"

        if not translating: # find start codon
            while True: 
                if i == len(self.seq): break # in case there is no start codon

                frame = frame[1:3] + self.seq[i]
                i += 1 

                if frame == "ATG": 
                    translating = True; protein = "M"
                    break

        if translating: 
            while True: 
                frame = self.seq[i:i+3]
                i += 3

                if i > len(self.seq): break

                aminoacid = CODONS[frame][1] # only want the 1-letter a.a. code
                if aminoacid == "O": break # do not translate past the first "stop" codon
                protein = protein + aminoacid
        return protein

if __name__ == "__main__":
    tseq1 = NucleicAcid("ATGGGCGC")
    tseq2 = NucleicAcid("AUGGGCGC", False) 
    tseq3 = NucleicAcid("ATCGatcg") 
    tseq4 = NucleicAcid("ATGCTGGTGACG")
    tseq5 = NucleicAcid("AAAAAAAATGCTGACGTGA")
    tseq6 = NucleicAcid("AAAAAAAAGTGA")
    tseq7 = NucleicAcid("ATGCTGGTGACGA")

    assert tseq1.rev_comp() == "GCGCCCAT" 
    assert tseq2.rev_comp() == "GCGCCCAU" 
    print("Reverse complement working correctly.")
    
    assert tseq3.gc_content() == 0.5
    assert tseq4.gc_content() == 7/12
    print("GC content correctly calculated.")

    tseq1.transcribe()
    assert tseq1.seq == "AUGGGCGC"
    assert tseq1.is_dna == False 
    tseq1.rev_transcribe()
    assert tseq1.seq == "ATGGGCGC"
    assert tseq1.is_dna == True 
    tseq4.transcribe()
    assert tseq4.seq == "AUGCUGGUGACG"
    assert tseq4.is_dna == False 
    tseq4.rev_transcribe()
    assert tseq4.seq == "ATGCTGGTGACG"
    assert tseq4.is_dna ==True 
    print("Transcription and reverse transcription working correctly.")

    assert tseq4.translate() == "MLVT"
    assert tseq5.translate() == "MLT"
    assert tseq6.translate() == ""
    assert tseq7.translate() == "MLVT"
    print("Translation working correctly.")
