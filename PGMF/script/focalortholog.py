#!/usr/bin/env python

""" 
Purpse:
    generate ortholog summary
    (pairwise between Dmel and non-melanogaster)
"""
from sharedinfo import get_A2B, get_M2N, DxxxGID_to_YOgnID, ordered_species
from focalannotation import FocalAnnotation

class FocalOrtholog:
    """FocalOrtholog object"""
    def __init__(self, species):
        self.species = species

        # ID conversion between annotation
        self.DxxxGID2FBgnID = get_M2N("../data/jaccard/" + species + ".jaccard", 0, 1)
        self.FBgnID2DxxxGID = get_M2N("../data/jaccard/" + species + ".jaccard", 1, 0)
        
        # known one2one ortholog concise
        self.FBgnID2dmelFBgnID = get_M2N("../data/ortholog/" + species + ".concise.one2one", 2, 0)

        # neo one2one ortholog
        self.DxxxGID2dmelFBgnID = get_M2N("../data/ortholog/" + species + ".neo121", 1, 0)
        self.dmelFBgnID2DxxxGID = get_M2N("../data/ortholog/" + species + ".neo121", 0, 1)

def generate_ortholog_table_for_species(species):
    o = FocalOrtholog(species)
    with open("../data/ortholog/" + species + ".ortholog.txt", "w") as f:
        f.write("\t".join(["YOgnID", "FBgnID", "Dmel"]) + "\n")
        for FBgnID in o.FBgnID2dmelFBgnID.keys():
            dmelFBgnID = o.FBgnID2dmelFBgnID[FBgnID]
            
            """ connecting by annotation-connection """
            try:
                DxxxGID = o.FBgnID2DxxxGID[FBgnID]
            except:
                f.write("-" + "\t" + FBgnID + "\t" + dmelFBgnID + "\n")
            else:
                f.write(DxxxGID_to_YOgnID(DxxxGID) + "\t" +  FBgnID + "\t" + dmelFBgnID + "\n")
                    
        for DxxxGID in o.DxxxGID2dmelFBgnID.keys():
            dmelFBgnID = o.DxxxGID2dmelFBgnID[DxxxGID]
            f.write(DxxxGID_to_YOgnID(DxxxGID) + "\t" + "-" + "\t" + dmelFBgnID + "\n")

def get_all_ortholog_tables():
    for species in ordered_species:
        generate_ortholog_table_for_species(species)

def join_all_ortholog_tables():
    dct = dict()
    for species in ordered_species:
        dct0 = get_M2N("../data/ortholog/" + species + ".ortholog.txt", 2, 0)
        dct1 = get_M2N("../data/ortholog/" + species + ".ortholog.txt", 2, 1)
        for dmelFBgnID in dct0.keys():
            YOgnID = dct0[dmelFBgnID]
            FBgnID = dct1[dmelFBgnID]
            if not dmelFBgnID in dct.keys():
                dct[dmelFBgnID] = dict()
            dct[dmelFBgnID][species] = [YOgnID, FBgnID]

    count = 0
    with open("../data/ortholog/ortholog.txt", "w") as f:
        f.write("\t".join(["FBgnID_Dmel", "GeneSymbol_Dmel"]) + "\t" + "\t".join(ordered_species) + "\n")
        for dmelFBgnID in dct.keys():
            count += 1
            print(count)
            dmelFBgnID2dmelSymbol = FocalAnnotation("dmel").dmelFBgn2dmelSymbol
            dmelSymbol = dmelFBgnID2dmelSymbol[dmelFBgnID]
            f.write(dmelFBgnID + "\t" + dmelSymbol)
            for species in ordered_species:
                ID = "-"
                if species in dct[dmelFBgnID].keys():
                    (YOgnID, FBgnID) = dct[dmelFBgnID][species]
                    if YOgnID == "-":
                        ID = FBgnID
                    else:
                        ID = YOgnID
                f.write("\t" + ID)
            f.write("\n")
                
            
if __name__ == '__main__':
    # get_all_ortholog_tables()
    # join_all_ortholog_tables()
