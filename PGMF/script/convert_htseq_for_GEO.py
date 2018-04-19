#!/usr/bin/env python

"""
Purpose: convert htseq matrix for GEO
"""
from __future__ import print_function, division
from sharedinfo import get_lines, DxxxGID_to_YOgnID, jaccard
from focalannotation import get_id


class Htseq:
    """ Htseq object """
    def __init__(self, sample):
        self.sample = sample
        self.species, self.sex, self.tissue, self.replicate = self.sample.split("_")
        self.v3exonmap = self.get_exonmap("v3")
        self.fbexonmap = self.get_exonmap("fb")
        self.v3htseqlines = get_lines("../data/htseq/FB2017_03_v3/", sample + ".htseq_reverse_HiSAT2.txt")
        self.fbhtseqlines = get_lines("../data/htseq/FB2017_03/", sample + ".htseq_reverse_HiSAT2.txt")
        self.DxxxGID2FBgnID = self.get_DxxxGID2FBgnID()
        self.DxxxGID2jaccard = self.get_DxxxGID2jaccard()
        #self.convert_DxxxGID_to_YOgnID_in_v3htseq()
        
    def convert_DxxxGID_to_YOgnID_in_htseq(self):
        nlines = list()
        for line in self.v3htseqlines:
            if not line.startswith("__"):
                DxxxGID, rc = line.rstrip().split("\t")
                nlines.append(DxxxGID_to_YOgnID(DxxxGID) + "\t" + rc + "\n")
        self.v3htseqlines = nlines
 
    def get_DxxxGID2FBgnID(self):
        """ get connection between DxxxGID and FBgn """
        lines = get_lines("../data/ortholog/", self.species + ".ee")
        DxxxGID2FBgnID = dict()
        for line in lines:
            temp = line.rstrip().split()
            DxxxGID = temp[0]
            FBgnID = temp[1]
            DxxxGID2FBgnID[DxxxGID] = FBgnID
        return(DxxxGID2FBgnID)

    def get_exonmap(self, tag):
        if tag == "v3":
            fn = "../data/annotation/" + self.species + ".v3.gtf"
        elif tag == "fb":
            fn = "../data/annotation/" + self.species + ".gtf"
        dct = dict()
        with open(fn, "r") as f:
            for line in f.readlines():
                (scaffold, tag, feature, start, end, scoare, strand, dot, others) = line.rstrip().split("\t")
                this_geneid = get_id(others, 'gene_id')
                if not this_geneid in dct.keys():
                    dct[this_geneid] = set()
                dct[this_geneid].add(start + "." + end)
        return(dct)

    def get_DxxxGID2jaccard(self):
        """ get jaccard between connection DxxxGID and FBgnID """
        DxxxGID2jaccard = dict()
        count = 0
        for DxxxGID in self.DxxxGID2FBgnID.keys():
            count += 1
            if count % 100 == 0:
                print(count)
            FBgnID = self.DxxxGID2FBgnID[DxxxGID]
            v3exonmap = self.v3exonmap[DxxxGID]
            fbexonmap = self.fbexonmap[FBgnID]
            j = jaccard(v3exonmap, fbexonmap)
            # print(DxxxGID, FBgnID, j)
            DxxxGID2jaccard[DxxxGID] = j
        return(DxxxGID2jaccard)

if __name__ == '__main__':
    h = Htseq("w1118_m_wb_R1")    
