#!/usr/bin/env python

""" 
Purpse:
    generate ortholog summary
    (pairwise between Dmel and non-melanogaster)
"""
from sharedinfo import get_A2B, get_M2N, DxxxGID_to_YOgnID

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

if __name__ == '__main__':
    o = FocalOrtholog("dyak")
    for FBgnID in o.FBgnID2dmelFBgnID.keys():
        dmelFBgnID = o.FBgnID2dmelFBgnID[FBgnID]
        
        """ connecting by annotation-connection """
        try:
            DxxxGID = o.FBgnID2DxxxGID[FBgnID]
        except:
            print("-" + "\t" + FBgnID + "\t" + dmelFBgnID)
        else:
            print(DxxxGID_to_YOgnID(DxxxGID) + "\t" +  FBgnID + "\t" + dmelFBgnID)
                
    for DxxxGID in o.DxxxGID2dmelFBgnID.keys():
        dmelFBgnID = o.DxxxGID2dmelFBgnID[DxxxGID]
        print(DxxxGID_to_YOgnID(DxxxGID) + "\t" + "-" + "\t" + dmelFBgnID)
