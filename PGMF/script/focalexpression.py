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
        
        # v2 updated annotation
        self.filename = species + ".expression.nrc.tab"
        self.lines = get_lines("../data/expression", self.filename)
        self.expressionheader, self.geneid2expression = self.get_geneid2expression(1)

        # FB2017_03 old annotation
        self.oldfilename = species + ".expression.onrc.tab"
        self.oldlines = get_lines("../data/expression", self.oldfilename)
        self.oldexpressionheader, self.oldgeneid2expression = self.get_geneid2expression(0)

    def get_geneid2expression(self, version):
        """
           get expression matrix 
           version = 0 or 1 (old or updated)
        """
        dct_geneid2expression = dict()

        if version == 0:
            lines = self.oldlines
        elif version == 1:
            lines = self.lines

        expressionheader = lines[0]
        for line in lines[1:]:
            elements = line.rstrip().split("\t")
            geneid, expression = elements[0], elements[1:]
            dct_geneid2expression[geneid] = expression
        return expressionheader, dct_geneid2expression
    

    
if __name__ == '__main__':
    # """ generate network expression matrix for yijie """
    # for dxxx in sharedinfo.ordered_species: 
    #    ann = focalannotation.FocalAnnotation(dxxx)
    #    exp = FocalExpression(dxxx)
    #    with open("../data/output/" + dxxx + ".genic.nrc.txt", 'w') as f:
    #        f.write(exp.expressionheader)
    #        for geneid in exp.geneid2expression.keys():
    #            expression = exp.geneid2expression[geneid]
    #            if geneid in ann.geneid2refgeneid121.keys():
    #                refgeneid = ann.geneid2refgeneid121[geneid]
    #                f.write(refgeneid + "\t" + "\t".join(expression) + "\n")

    """ generate network expression matrix for yijie (orthologs converted)"""
    dct_ortholog_expression = dict()
    for dxxx in sharedinfo.ordered_species:
        print ("Now processing " + dxxx)
        ann = focalannotation.FocalAnnotation(dxxx)
        exp = FocalExpression(dxxx)
        with open("../data/output/" + dxxx + ".genic.nrc.txt", 'w') as f:
            f.write(exp.expressionheader)
            for geneid in exp.geneid2expression.keys():
                expression = exp.geneid2expression[geneid]
                if geneid in ann.geneid2refgeneid121.keys():
                    refgeneid = ann.geneid2refgeneid121[geneid]
                    #f.write(refgeneid + "\t" + "\t".join(expression) + "\n")

                    findDmelOrtholog = 0
                    if refgeneid.endswith(dxxx):
                        dmelrefgeneid = refgeneid.replace("_" + dxxx, "")
                        findDmelOrtholog += 1
                    elif refgeneid in ann.olo121.keys():
                        dmelrefgeneid = ann.olo121[refgeneid]
                        findDmelOrtholog += 1

                    if findDmelOrtholog > 0:
                        if not dmelrefgeneid in dct_ortholog_expression.keys():
                            dct_ortholog_expression[dmelrefgeneid] = dict()
                        f.write(dmelrefgeneid + "\t" + "\t".join(expression) + "\n")
                        dct_ortholog_expression[dmelrefgeneid][dxxx] = expression

    #with open("../data/output/dmel.ortholog.expression.longlist", 'w') as f:
    #    for dmelrefgeneid in dct_ortholog_expression.keys():
    #        for dxxx in dct_ortholog_expression[dmelrefgeneid]:
    #            f.write(dmelrefgeneid + "\t" + dxxx + "\t" + "\t".join(dct_ortholog_expression[dmelrefgeneid][dxxx]) + "\n")



    """ gene expression matrix based on old annotation """
    for dxxx in sharedinfo.ordered_species:
        print ("Now processing " + dxxx)
        ann = focalannotation.FocalAnnotation(dxxx)
        exp = FocalExpression(dxxx)
        with open("../data/output/" + dxxx + ".genic.onrc.txt", 'w') as f:
            f.write(exp.oldexpressionheader)
            for geneid in exp.oldgeneid2expression.keys():
                expression = exp.oldgeneid2expression[geneid]
                if geneid in ann.olo121.keys():
                    dmelgeneid = ann.olo121[geneid]
                    f.write(dmelgeneid + "\t" + "\t".join(expression) + "\n")

    



    """ test """
    #for geneid in ann.geneid2refgeneid121.keys():
    #    if ann.geneid2refgeneid121[geneid] == "FBgn0000504_dyak":
    #        print(geneid, ann.geneid2refgeneid121[geneid])

    #MFBST.8049 FBgn0000504_dyak
    #MFBST.8207 FBgn0000504_dyak
