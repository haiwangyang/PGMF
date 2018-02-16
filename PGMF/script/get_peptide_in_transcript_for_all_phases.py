#!/usr/bin/env python

"""
   get all possible peptide in transcript
   the peptides include all three phases, and stop codons
"""

import re
import focalgene
from pyfaidx import Fasta

def get_indices_of_matches(seq, match):
    """ string: "Thisxis", match: "is", return [2, 5] """
    if match == "*":
        match = "\\" + match
    iter = re.finditer(match, seq)
    indices = [m.start(0) for m in iter]
    return(indices)

def get_indices_located_within_range(start, end, indices):
    """ start: 3, end: 10, indices: [1, 5, 9, 12], return [5, 9] """
    within = list()
    for i in range(start, end + 1):
        if i in indices:
            within.append(i)
    return within


def get_all_possible_protein(seq):
    """ seq is a long protein sequence (with many M and *)
        the function is trying to find all short peptide within it that make sense
        for instance,
        seq: SMSSSMSS*S
        return: {1-8:MSSSMSS*, 5-8:MSS*]
    """
    indices_M = get_indices_of_matches(seq, "M")
    indices_stop = get_indices_of_matches(seq, "*")
    start = 0
    dct = dict()
    for i_s in indices_stop:
        indices_M_within = get_indices_located_within_range(start, i_s, indices_M)
        #print(seq)
        #print(indices_M)
        #print(indices_stop)
        #print("M within " + str(start) + " and " + str(i_s))
        #print(indices_M_within)
        for i_M in indices_M_within:
            peptide = seq[i_M:i_s + 1]
            k = str(i_M) + "-" + str(i_s) 
            dct[k] = peptide
        #print("\n\n")
        start = i_s    
    return(dct)

if __name__ == '__main__':
    fasta = Fasta("../data/pacbio/pacbio_new_gene_model.fasta", duplicate_action="longest")
    with open("../data/pacbio/pacbio_new_gene_model.all_phase_peptide.fasta", "w") as f:
        for name in fasta.keys():
            seq = str(fasta[name])
            for i in range(3):
                dct_all_peptide = get_all_possible_protein(focalgene.get_pep_and_leftover_from_dna(seq, i)[0])
                for r in dct_all_peptide.keys():
                    f.write(">" + name + "=" + str(i) + "=" + r + "\n" + dct_all_peptide[r] + "\n")
