import re
from pyfaidx import Fasta

def main():
    """ read fasta of 26bp seq of 3' pacbio transcript """
    fasta = Fasta("../data/pacbio/pacbio_new_gene_model.bam.down26.fasta", duplicate_action="longest")
    for name in fasta.keys():
        seq = str(fasta[name])
        m = re.search('^(' + "A" + '+)', seq)
        if m:
            p = str(m.group(1))
            print(p + "\t" + str(len(p)))
        else:
            print("N" + "\t" + str(0))

if __name__ == '__main__':
    main()
