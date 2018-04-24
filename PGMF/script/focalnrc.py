#!/usr/bin/env python

"""
Purpose: convert DESeq2 normalized read count matrix for Supplementary Data
"""
from __future__ import print_function, division
from sharedinfo import get_lines, DxxxGID_to_YOgnID, jaccard, typical_strains, species2dxxx, sample2GSM, samples, ordered_species
from focalannotation import get_id
import os
import pandas


def convert_with_dict(lst, dct):
    """ convert lst to another lst
        with dct table if possible """
    lst2 = []
    for i in lst:
        try:
            ii = dct[i]
        except:
            ii = "-"
        else:
            pass
        lst2.append(ii)
    return(lst2)
                

class Nrc:
    """ Nrc object """
    def __init__(self, species):        
        self.species = species
        self.v3exonmap = self.get_exonmap("v3")
        self.fbexonmap = self.get_exonmap("fb")
        self.v3_nrc_table = pandas.read_table("../data/nrc/FB2017_03_v3/" + species + ".expression.nrc.txt")
        self.fb_nrc_table = pandas.read_table("../data/nrc/FB2017_03/" + species + ".expression.nrc.txt")
        #self.fbhtseq = self.get_fbhtseq()
        self.DxxxGID2YOgnID = self.get_DxxxGID2YOgnID()
        self.DxxxGID2FBgnID = self.get_DxxxGID2FBgnID()
        self.DxxxGID2jaccard = self.get_DxxxGID2jaccard()
        self.update_v3_nrc_table()

    def update_v3_nrc_table(self):
        """ update htseq txt """
        DxxxGIDs = list(self.v3_nrc_table.index)
        YOgnIDs = convert_with_dict(DxxxGIDs, self.DxxxGID2YOgnID)
        FBgnIDs = convert_with_dict(DxxxGIDs, self.DxxxGID2FBgnID)
        js = convert_with_dict(DxxxGIDs, self.DxxxGID2jaccard)
        self.v3_nrc_table.insert(loc = 0, column = "FBgnID", value = FBgnIDs)
        self.v3_nrc_table.insert(loc = 0, column = "jaccard", value = js)
        self.v3_nrc_table.insert(loc = 0, column = "YOgnID", value = YOgnIDs)
        self.v3_nrc_table.index = self.v3_nrc_table.YOgnID
        self.v3_nrc_table = self.v3_nrc_table.drop(["YOgnID"], axis=1)

    def print_out_updated_v3_nrc_table(self):
        if self.species == "dgri":
            self.v3_nrc_table[list(self.v3_nrc_table.columns.values)[0:66]].to_csv("../data/nrc/" + self.species + ".nrc.txt", sep="\t")
        else:        
            self.v3_nrc_table.to_csv("../data/nrc/" + self.species + ".nrc.txt", sep="\t")

    def get_DxxxGID2YOgnID(self):
        nlines = list()
        DxxxGID2YOgnID = dict()
        for line in list(self.v3_nrc_table.index):
            if not line.startswith("__"):
                temp = line.rstrip().split("\t")
                DxxxGID = temp[0]
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

def generate_new_nrc_for_all_species():
    for species in ordered_species:
        print("generating new nrc table for " + species)
        n = Nrc(species)
        n.print_out_updated_v3_nrc_table()

def main():
    # generate_new_nrc_for_all_species()
    n = Nrc("dgri")
    
if __name__ == '__main__':
    # main()
    n = Nrc("dgri")
