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
codons = {'AAA': ['Lys', 'K',	'Lysine'],
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

class nuc_acid(): 
    '''Nucleic acid sequence, default DNA, may be RNA.
    If RNA, user must define. 
    No 'N' bases allowed.''' 

    DEFAULT_KIND = "DNA"

    def __init__(self, sequence: str, identity: str = DEFAULT_KIND):
        self.kind: str = identity
        self.seq: str = sequence.upper()
        self.length: int = len(sequence)
        assert self.validateseq()

    def validateseq(self): 
        '''Must be a sequence containing only A T C or G if DNA,
        or a A U C or G if RNA.
        Case insensitive.'''
        if self.kind == "DNA": return set(DNA_bases).issuperset(self.seq)
        if self.kind == "RNA": return set(RNA_bases).issuperset(self.seq)

    def gc_content(self) -> float:
        '''Calculate GC content of sequence.'''
        gcount: int = self.seq.count('G') + self.seq.count('g')
        ccount: int = self.seq.count('C') + self.seq.count('c')
        return (ccount + gcount) / self.length

    def rev_comp(self) -> str:
        '''Returns the reverse complement of seq.'''
        reversecomp: str = ""
        if self.kind == "DNA":
            for index in range(len(self.seq)):
                reversecomp = reversecomp + DNA_comp[self.seq[-1*(index + 1)]]
        if self.kind == "RNA": 
            for index in range(len(self.seq)):
                reversecomp = reversecomp + RNA_comp[self.seq[-1*(index + 1)]]
        return reversecomp

    def transcribe(self):
        '''Transcribes DNA sequence into RNA sequence in-place. 
        Assumes DNA sequence is coding.
        This method changes seq to RNA and kind to "RNA"'''

        assert self.kind == "DNA"
        rna = ""
        for i in range(len(self.seq)): 
            if self.seq[i] == "T": rna = rna + "U"
            else: rna = rna + self.seq[i]
        self.seq = rna
        self.kind = "RNA"

    def rev_transcribe(self):
        '''Reverse-transcribes RNA sequence into DNA sequence in-place. 
        Seq becomes coding DNA sequence, type becomes "DNA"'''

        assert self.kind == "RNA"
        dna = ""
        for i in range(len(self.seq)): 
            if self.seq[i] == "U": dna = dna + "T"
            else: dna = dna + self.seq[i]
        self.seq = dna
        self.kind = "DNA"

    def translate(self) -> str: 
        '''Translate DNA nuc_acid. 
        Sequence must have a start codon to find the reading frame.
        Translates up to first stop codon encountered in seq.'''

        assert self.kind == "DNA"
        
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

                aminoacid = codons[frame][1] # only want the 1-letter a.a. code
                if aminoacid == "O": break # do not translate past the first "stop" codon
                protein = protein + aminoacid
        return protein

if __name__ == "__main__":
    tseq1 = nuc_acid("ATGGGCGC")
    tseq2 = nuc_acid("AUGGGCGC", "RNA") 
    tseq3 = nuc_acid("ATCGatcg") 
    tseq4 = nuc_acid("ATGCTGGTGACG")
    tseq5 = nuc_acid("AAAAAAAATGCTGACGTGA")
    tseq6 = nuc_acid("AAAAAAAAGTGA")
    tseq7 = nuc_acid("ATGCTGGTGACGA")

    assert tseq1.rev_comp() == "GCGCCCAT" 
    assert tseq2.rev_comp() == "GCGCCCAU" 
    print("Reverse complement working correctly.")
    
    assert tseq3.gc_content() == 0.5
    assert tseq4.gc_content() == 7/12
    print("GC content correctly calculated.")

    tseq1.transcribe()
    assert tseq1.seq == "AUGGGCGC"
    assert tseq1.kind == "RNA"
    tseq1.rev_transcribe()
    assert tseq1.seq == "ATGGGCGC"
    assert tseq1.kind == "DNA"
    tseq4.transcribe()
    assert tseq4.seq == "AUGCUGGUGACG"
    assert tseq4.kind == "RNA"
    tseq4.rev_transcribe()
    assert tseq4.seq == "ATGCTGGTGACG"
    assert tseq4.kind == "DNA"
    print("Transcription and reverse transcription working correctly.")

    assert tseq4.translate() == "MLVT"
    assert tseq5.translate() == "MLT"
    assert tseq6.translate() == ""
    assert tseq7.translate() == "MLVT"
    print("Translation working correctly.")
