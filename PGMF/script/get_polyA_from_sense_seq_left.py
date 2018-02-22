from focalisoseq import seq2polyA
from pyfaidx import Fasta

def main():
    """ read fasta of 26bp seq of 3' pacbio transcript """
    fasta = Fasta("../data/pacbio/pacbio_new_gene_model.bam.down26.fasta", duplicate_action="longest")
    for name in fasta.keys():
        seq = fasta[name]
        print(name)
        print(seq2polyA(seq))

if __name__ == '__main__':
    main()
