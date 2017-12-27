#!/usr/bin/env python

"""
Purpose:
    Handling normalized read count of genes

"""
import focalannotation
import sharedinfo
from sharedinfo import exist_file, get_lines

class FocalExpression:
    """FocalExpression object"""
    def __init__(self, species):
        self.species = species
        self.filename = species + ".expression.nrc.tab"
        self.lines = get_lines("../data/expression", self.filename)
        self.expressionheader, self.geneid2expression = self.get_geneid2expression()

    def get_geneid2expression(self):
        """ get expression matrix """
        dct_geneid2expression = dict()
        expressionheader = self.lines[0]
        for line in self.lines[1:]:
            elements = line.rstrip().split("\t")
            geneid, expression = elements[0], elements[1:]
            dct_geneid2expression[geneid] = expression
        return expressionheader, dct_geneid2expression
    

    
if __name__ == '__main__':
    """ generate network expression matrix for yijie """
    for dxxx in sharedinfo.ordered_species: 
        ann = focalannotation.FocalAnnotation(dxxx)
        exp = FocalExpression(dxxx)
        with open("../data/output/" + dxxx + ".genic.nrc.txt", 'w') as f:
            f.write(exp.expressionheader)
            for geneid in exp.geneid2expression.keys():
                expression = exp.geneid2expression[geneid]
                if geneid in ann.geneid2refgeneid121.keys():
                    refgeneid = ann.geneid2refgeneid121[geneid]
                    f.write(refgeneid + "\t" + "\t".join(expression) + "\n")

    """ test """
    #for geneid in ann.geneid2refgeneid121.keys():
    #    if ann.geneid2refgeneid121[geneid] == "FBgn0000504_dyak":
    #        print(geneid, ann.geneid2refgeneid121[geneid])

    #MFBST.8049 FBgn0000504_dyak
    #MFBST.8207 FBgn0000504_dyak