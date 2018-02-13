#!/usr/bin/env python

"""
   get all possible peptide in transcript
   the peptides include all three phases, and stop codons
"""

import focalgene
from pyfaidx import Fasta

if __name__ == '__main__':
    fasta = Fasta("../data/pacbio/pacbio_new_gene_model.fasta", duplicate_action="longest")
    with open("../data/pacbio/pacbio_new_gene_model.all_phase_peptide.fasta", "w") as f:
        for name in fasta.keys():
            seq = str(fasta[name])
            for i in range(3):
                f.write(">" + name + "=" + str(i) + "\n" + focalgene.get_pep_and_leftover_from_dna(seq, i)[0] + "\n")



