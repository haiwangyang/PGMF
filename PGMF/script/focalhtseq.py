#!/usr/bin/env python

"""
Purpose: convert htseq matrix for GEO
"""
from __future__ import print_function, division
from sharedinfo import get_lines, DxxxGID_to_YOgnID, jaccard, typical_strains, species2dxxx, sample2GSM, samples
from focalannotation import get_id
import os

class Htseq:
    """ Htseq object """
    def __init__(self, sample):
        self.sample = sample
        self.species, self.sex, self.tissue, self.replicate = self.sample.split("_")
        self.v3exonmap = self.get_exonmap("v3")
        self.fbexonmap = self.get_exonmap("fb")
        self.v3htseqlines = get_lines("../data/htseq/FB2017_03_v3/", sample + ".htseq_reverse_HiSAT2.txt")
        self.fbhtseqlines = get_lines("../data/htseq/FB2017_03/", sample + ".htseq_reverse_HiSAT2.txt")
        self.fbhtseq = self.get_fbhtseq()
        self.DxxxGID2YOgnID = self.get_DxxxGID2YOgnID()
        self.DxxxGID2FBgnID = self.get_DxxxGID2FBgnID()
        self.DxxxGID2jaccard = self.get_DxxxGID2jaccard()
        self.update_v3htseq_lines()


    def update_v3htseq_lines(self):
        """ update htseq txt """
        GSMinfos = sample2GSM[self.sample]
        for GSMinfo in GSMinfos:
            GSM, sampleInGEO = GSMinfo
            with open("../data/htseq/" + GSM + "_" + sampleInGEO + ".htseq_reverse.txt", "w") as f:
                for line in self.v3htseqlines:
                    if not line.startswith("__"):
                        DxxxGID, rc = line.rstrip().split("\t")
                        YOgnID = self.DxxxGID2YOgnID[DxxxGID]
                        try:
                            FBgnID = self.DxxxGID2FBgnID[DxxxGID]
                        except:
                            FBgnID = "-"
                            j = "-"
                            fbexp = "-"
                        else:
                            FBgnID = self.DxxxGID2FBgnID[DxxxGID]
                            j = self.DxxxGID2jaccard[DxxxGID]
                            fbexp = self.fbhtseq[FBgnID]
                    f.write("\t".join([YOgnID, str(j), FBgnID, rc]) + "\n")

    def get_fbhtseq(self):
        """ get fb htseq in dct """
        dct = dict()
        for line in self.fbhtseqlines:
            if not line.startswith("__"):
                FBgnID, rc = line.rstrip().split("\t")
                dct[FBgnID] = rc
        return(dct)

    def get_DxxxGID2YOgnID(self):
        nlines = list()
        DxxxGID2YOgnID = dict()
        for line in self.v3htseqlines:
            if not line.startswith("__"):
                DxxxGID, rc = line.rstrip().split("\t")
                DxxxGID2YOgnID[DxxxGID] = DxxxGID_to_YOgnID(DxxxGID)
        return(DxxxGID2YOgnID)

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
        dxxx = species2dxxx[self.species]
        jfile = "../data/jaccard/" + dxxx + ".jaccard"
        try:
            os.stat(jfile)
        except:
            count = 0
            for DxxxGID in self.DxxxGID2FBgnID.keys():
                count += 1
                if count % 1000 == 0:
                    print(count)
                FBgnID = self.DxxxGID2FBgnID[DxxxGID]
                v3exonmap = self.v3exonmap[DxxxGID]
                fbexonmap = self.fbexonmap[FBgnID]
                j = jaccard(v3exonmap, fbexonmap)
                # print(DxxxGID, FBgnID, j)
                DxxxGID2jaccard[DxxxGID] = j
        else:
            with open(jfile, "r") as f:
                for line in f.readlines():
                    DxxxGID, FBgnID, j = line.rstrip().split("\t")
                    DxxxGID2jaccard[DxxxGID] = j 
        return(DxxxGID2jaccard)

    

    def generate_jaccard_file(self):
        species = self.species
        if self.species == "w1118" or self.species == "oreR":
            species = "dmel"
        if self.species == "dgriG1":
            species = "dgri"
        with open("../data/jaccard/" + species + ".jaccard", "w") as f:
            for DxxxGID in self.DxxxGID2jaccard.keys():
                FBgnID = self.DxxxGID2FBgnID[DxxxGID]
                j = self.DxxxGID2jaccard[DxxxGID]
                f.write("\t".join([DxxxGID, FBgnID, str(j)]) + "\n")

def generate_jaccard_file_all_species():
    for s in typical_strains:
        print(s)
        h = Htseq(s + "_m_wb_R1")
        h.generate_jaccard_file()

def generate_GEO_htseq_txt_for_all_samples():
    for sample in samples:
        print(sample)
        h = Htseq(sample)

def main():
    # generate_jaccard_file_all_species()
    # h = Htseq("w1118_m_wb_R1")
    generate_GEO_htseq_txt_for_all_samples()

if __name__ == '__main__':
    main()    
