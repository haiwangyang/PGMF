#!/usr/bin/env python

"""
Purpose:
    Handling signalp report
    https://hpc.nih.gov/apps/signalp.html
"""


from __future__ import print_function, division
from sharedinfo import exist_file, get_lines
from pyfaidx import Fasta

def get_position2manualinfo(filename):
    """ position is the basic concised position from isoseqid
        3L:15453948-15455098-
        manuall info: sample, intron_containing, coding, conservation, intergenic
    """
    dct = dict()
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            (sample, intron, isoseqid, coding, conservation, intergenic) = line.rstrip().split("\t")
            position = isoseqid.split(".")[0]
            dct[position] = [sample, intron, isoseqid, coding, conservation, intergenic]
    return(dct)

def positionplus2position(positionplus):
    """ positionplus is 3R:863419-864511(-)=0=332-341
        position is 3R:863420-864511- (bed format)
    """
    positionraw = positionplus.split("=")[0]
    chrom, other = positionraw.split(":")
    start_end, strand = other.split("(")
    strand = strand.replace(")", "")
    start, end = start_end.split("-")
    return(chrom + ":" + str(int(start) + 1) + "-" + end + strand)

class FocalSignalp:
    """FocalSignalp object"""
    def __init__(self, filename):
        """ filename example: pacbio_new_gene_model.all_phase_peptide """
        self.filename = filename
        self.peptidefasta = Fasta("../data/pacbio/" + filename + ".fasta")
        print(self.peptidefasta)
        self.signalplines = get_lines("../data/pacbio", filename + ".fasta.signalp")
        self.position2manualinfo = get_position2manualinfo("../data/pacbio/pacbio_new_gene_model.tab")

    def print_out_characteristics_for_signalp(self):
        """ print out peptide length distribution, cleavage site, D value
            for peptide with signalp and without
        """
        with open("../data/pacbio/" + self.filename + ".fasta.signalp.summary", "w") as f:
            for line in self.signalplines:
                temp = line.rstrip().split("\t")
                if line.startswith("Name="):
                    positionplus = temp[0].replace("Name=", "")
                    position = positionplus2position(positionplus)
                    (sample, intron, isoseqid, coding, conservation, intergenic) = self.position2manualinfo[position]
                    peptide = str(self.peptidefasta[positionplus])
                    len_peptide = len(peptide)
                    if "SP='YES'" in line:
                        cleavage_site = temp[1].split(" ")[8]
                        D = temp[1].split(" ")[9].replace("D=", "")
                        f.write("\t".join(["signalp", positionplus, sample, intron, isoseqid, coding, conservation, intergenic, cleavage_site, str(len_peptide), D, peptide]) + "\n")
                    elif "SP='NO'" in line:
                        D = "NA"
                        cleavage_site = "NA"
                        f.write("\t".join(["not_signalp", positionplus, sample, intron, isoseqid, coding, conservation, intergenic, cleavage_site, str(len_peptide), D, peptide]) + "\n")
    
 
if __name__ == '__main__':
    fs = FocalSignalp("pacbio_new_gene_model.all_phase_peptide")
    fs.print_out_characteristics_for_signalp()    
