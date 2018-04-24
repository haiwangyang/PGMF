#!/usr/bin/env python

"""
Purpose: convert htseq matrix for GEO
"""
from __future__ import print_function, division
from sharedinfo import get_lines, DxxxTID_to_YOtrID, jaccard, typical_strains, species2dxxx, sample2GSM, samples
from focalannotation import get_id
import os

class Salmon:
    """ Salmon object """
    def __init__(self, sample):
        self.sample = sample
        self.species, self.sex, self.tissue, self.replicate = self.sample.split("_")
        self.v3exonmap = self.get_exonmap("v3") # transcript level
        self.fbexonmap = self.get_exonmap("fb") # transcript level
        self.v3salmonlines = get_lines("../data/salmon/FB2017_03_v3/", sample + ".salmon.txt")
        self.fbsalmonlines = get_lines("../data/salmon/FB2017_03/", sample + ".salmon.txt")
        self.fbsalmon = self.get_fbsalmon()
        self.DxxxTID2YOtrID = self.get_DxxxTID2YOtrID()
        #self.DxxxTID2FBtrID = self.get_DxxxTID2FBtrID()
        #self.DxxxTID2jaccard = self.get_DxxxTID2jaccard()
        #self.update_v3salmon_lines()

    def update_v3salmon_lines(self):
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

    def get_fbsalmon(self):
        """ get fb salmon in dct """
        dct = dict()
        for line in self.fbsalmonlines[1:]:
            if not line.startswith("__"):
                temp = line.rstrip().split("\t")
                dct[temp[0]] = temp[3]
        return(dct)

    def get_DxxxTID2YOtrID(self):
        nlines = list()
        DxxxTID2YOtrID = dict()
        for line in self.v3salmonlines[1:]:
            if not line.startswith("__"):
                temp = line.rstrip().split("\t")
                DxxxTID2YOtrID[temp[0]] = DxxxTID_to_YOtrID(temp[0])
        return(DxxxTID2YOtrID)

    def get_DxxxTID2FBtrID(self):
        """ get connection between DxxxTID and FBtr """
        lines = get_lines("../data/ortholog/", self.species + ".tee")
        DxxxGID2FBgnID = dict()
        for line in lines:
            temp = line.rstrip().split()
            DxxxTID = temp[1]
            FBtrID = temp[3]
            DxxxTID2FBtrID[DxxxTID] = FBtrID
        return(DxxxTID2FBtrID)

    def get_exonmap(self, tag):
        if tag == "v3":
            fn = "../data/annotation/" + self.species + ".v3.gtf"
        elif tag == "fb":
            fn = "../data/annotation/" + self.species + ".gtf"
        dct = dict()
        with open(fn, "r") as f:
            for line in f.readlines():
                (scaffold, tag, feature, start, end, scoare, strand, dot, others) = line.rstrip().split("\t")
                this_transid = get_id(others, 'transcript_id')
                if not this_transid in dct.keys():
                    dct[this_transid] = set()
                dct[this_transid].add(start + "." + end)
        return(dct)

    def get_DxxxTID2jaccard(self):
        """ get jaccard between connection DxxxGID and FBgnID """
        DxxxTID2jaccard = dict()
        dxxx = species2dxxx[self.species]
        jfile = "../data/jaccard/" + dxxx + ".t.jaccard"
        try:
            os.stat(jfile)
        except:
            count = 0
            for DxxxGID in self.DxxxTID2FBtrID.keys():
                count += 1
                if count % 1000 == 0:
                    print(count)
                FBgnID = self.DxxxTID2FBtrID[DxxxTID]
                v3exonmap = self.v3exonmap[DxxxTID]
                fbexonmap = self.fbexonmap[FBtrID]
                j = jaccard(v3exonmap, fbexonmap)
                # print(DxxxTID, FBtrID, j)
                DxxxTID2jaccard[DxxxTID] = j
        else:
            with open(jfile, "r") as f:
                for line in f.readlines():
                    DxxxTID, FBtrID, j = line.rstrip().split("\t")
                    DxxxTID2jaccard[DxxxTID] = j 
        return(DxxxTID2jaccard)

    

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
    # for s in typical_strains:
    for s in ["dyak",]:
        print(s)
        h = Htseq(s + "_m_wb_R1")
        h.generate_jaccard_file()

def generate_GEO_htseq_txt_for_all_samples():
    for sample in samples:
        print(sample)
        h = Htseq(sample)

def main():
    generate_jaccard_file_all_species()
    # h = Htseq("w1118_m_wb_R1")
    # generate_GEO_htseq_txt_for_all_samples()

if __name__ == '__main__':
    # main()    
    s = Salmon("dyak_m_wb_R1")
