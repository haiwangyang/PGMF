#!/usr/bin/env python

"""
Purpose: convert htseq matrix for GEO
"""
from __future__ import print_function, division
from sharedinfo import get_lines, DxxxTID_to_YOtrID, jaccard, typical_strains, species2dxxx, sample2GSM, samples, ordered_species
from focalannotation import get_id
import os
import pandas

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
        self.DxxxTID2FBtrID = self.get_DxxxTID2FBtrID()
        self.DxxxTID2jaccard = self.get_DxxxTID2jaccard()
        # self.update_v3salmon_lines()

    def update_v3salmon_lines(self):
        """ update htseq txt """
        GSMinfos = sample2GSM[self.sample]
        for GSMinfo in GSMinfos:
            GSM, sampleInGEO = GSMinfo
            with open("../data/salmon/" + sampleInGEO + ".salmon.txt", "w") as f:
                for line in self.v3salmonlines[1:]:
                    if not line.startswith("__"):
                        temp = line.rstrip().split("\t")
                        DxxxTID = temp[0]
                        YOtrID = self.DxxxTID2YOtrID[DxxxTID]
                        try:
                            FBtrID = self.DxxxTID2FBtrID[DxxxTID]
                        except:
                            FBtrID = "-"
                            j = "-"
                            # fbexp = "-"
                        else:
                            FBtrID = self.DxxxTID2FBtrID[DxxxTID]
                            j = self.DxxxTID2jaccard[DxxxTID]
                            # fbexp = self.fbsalmon[FBtrID]
                        f.write("\t".join([YOtrID, str(j), FBtrID, temp[3]]) + "\n")

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
        DxxxTID2FBtrID = dict()
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
            for DxxxTID in self.DxxxTID2FBtrID.keys():
                count += 1
                if count % 1000 == 0:
                    print(count)
                FBtrID = self.DxxxTID2FBtrID[DxxxTID]
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
        with open("../data/jaccard/" + species + ".t.jaccard", "w") as f:
            for DxxxTID in self.DxxxTID2jaccard.keys():
                FBtrID = self.DxxxTID2FBtrID[DxxxTID]
                j = self.DxxxTID2jaccard[DxxxTID]
                f.write("\t".join([DxxxTID, FBtrID, str(j)]) + "\n")

def generate_jaccard_file_all_species():
    for s in typical_strains:
    #for s in ["dyak",]:
        print(s)
        h = Salmon(s + "_m_wb_R1")
        h.generate_jaccard_file()

def generate_GEO_salmon_txt_for_all_samples():
    for sample in samples:
        print(sample)
        s = Salmon(sample)
        s.update_v3salmon_lines()

def follow_nrc_matrix_format():
    for species in ordered_species:
        print("working on " + species)
        with open("../data/nrc/" + species + ".nrc.txt") as f:
            lines = f.readlines()
            samples = lines[0].rstrip().split("\t")[3:]
            sample0 = samples[0]
            pdt0 = pandas.read_table("../data/salmon/" + sample0 + ".salmon.txt", names = ["YOtrID", "jaccard", "FBtrID", "TPM"])            
            merged = pdt0[["YOtrID", "jaccard", "FBtrID"]]
            for sample in samples:
                pdt = pandas.read_table("../data/salmon/" + sample + ".salmon.txt", names = ["YOtrID", "jaccard", "FBtrID", "TPM"])
                extraone = pdt[["TPM"]].rename(columns={"TPM":sample})
                merged = merged.join(extraone)
            merged.to_csv("../data/salmon/" + species + ".transTPM.txt", sep="\t", index = False)

def main():
    # generate_jaccard_file_all_species()
    # s = Salmon("w1118_m_wb_R1")
    # generate_GEO_salmon_txt_for_all_samples()
    follow_nrc_matrix_format()    

if __name__ == '__main__':
    main()


