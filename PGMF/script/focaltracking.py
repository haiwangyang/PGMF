#!/usr/bin/env python

"""
Purpose:
    Handling tracking

"""
import sharedinfo
import focalisoseq
import re
from sharedinfo import exist_file, get_lines
from focalannotation import get_id

def print_overlaptype_intronnum(sample):
    species, sex, tissue, replicate = sample.split("_")
    ft = FocalTracking(species, sex, tissue, replicate)
    fi = focalisoseq.FocalIsoseq(species, sex, tissue, replicate)
    with open("../data/pacbio/jaccard_zero." + sample + ".exon.gtf.tracking2", "w") as f:
        for isoseqid in ft.isoseqid2overlaptype.keys():
            overlaptype = ft.isoseqid2overlaptype[isoseqid]
            position = fi.isoseqid2position[isoseqid]
            if isoseqid in fi.isoseqid2intronnum:
                intronnum = fi.isoseqid2intronnum[isoseqid]
                f.write(position + "\t" + isoseqid + "\t" + overlaptype + "\t" + str(intronnum) + "\n")

class FocalTracking:
    """FocalTracking object"""
    def __init__(self, species, sex, tissue, replicate):
        self.species = species
        self.sex = sex
        self.tissue = tissue
        self.replicate = replicate
        self.sample = species + "_" + sex + "_" + tissue + "_" + replicate
        self.geneidtransid2strand = self.get_geneidtransid2strand()
        self.trackinglines = get_lines("../data/pacbio/", "jaccard_zero." +  self.sample + ".exon.gtf.tracking")
        self.isoseqid2overlaptype = self.get_isoseqid2overlaptype()    

    def get_geneidtransid2strand(self):
        """ get transid 2 strand """
        dct = dict()
        for line in get_lines("../data/annotation/", self.species + ".v3.gtf"):
            (scaffold, tag, feature, start, end, scoare, strand, dot, others) = line.rstrip().split("\t")
            this_geneid = get_id(others, 'gene_id')
            this_transid = get_id(others, 'transcript_id')
            dct[this_geneid + "|" + this_transid] = strand
        return dct

    def get_isoseqid2overlaptype(self):
        """ get new gene from pacbio isoseq """
        dct = dict()
        for line in self.trackinglines:
            TCONS, XLOC, geneidtransid, type, other = line.rstrip().split()
            strand_geneidtransid = ""
            if geneidtransid in self.geneidtransid2strand.keys():
                strand_geneidtransid = self.geneidtransid2strand[geneidtransid]
            position_isoseqid = other.split("|")[1]
            position, isoseqid = position_isoseqid.split(".")

            strand_isoseqid = position[-1]

            if type == "i": # inside annotated intron
                #print (strand_geneidtransid, strand_isoseqid, line)
                if strand_geneidtransid != strand_isoseqid: # anti_strand inside annotated intron:
                    dct[isoseqid] = "pass_ia"
                else:                                       # sense_strand inside annotated intron:
                    dct[isoseqid] = "fail_is"             
            elif type == "p":
                dct[isoseqid] = "pass_p"
            elif type == "u":
                dct[isoseqid] = "pass_u"                
            else:
                dct[isoseqid] = "fail_" + type
        return dct


if __name__ == '__main__':
    for sample in sharedinfo.pacbio_sample:
        print_overlaptype_intronnum(sample)
