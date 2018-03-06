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
        
        # v3 updated annotation
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
    

def main():
    """ generate network expression matrix for yijie (orthologs converted)"""
    dmelFBgn2dmelSymbol = focalannotation.FocalAnnotation("dmel").dmelFBgn2dmelSymbol

    """ each species expression """
    dct_ortholog_expression = dict()
    dct_ortholog_header = dict()
    for dxxx in sharedinfo.ordered_species:
        print ("Now processing " + dxxx)
        ann = focalannotation.FocalAnnotation(dxxx)
        exp = FocalExpression(dxxx)
        with open("../data/expression/" + dxxx + ".genic.nrc.txt", 'w') as f:
            f.write("GeneID" + "\t" + "ortholog_dmelFBgn" + "\t" + "ortholog_dmelSymbol" + "\t" + exp.expressionheader)
            dct_ortholog_header[dxxx] = exp.expressionheader
            for GdxxxID in exp.geneid2expression.keys():
                expression = exp.geneid2expression[GdxxxID]
                if GdxxxID in ann.GdxxxID2dmelFBgn.keys():    
                    dmelFBgn = ann.GdxxxID2dmelFBgn[GdxxxID]
                    dmelSymbol = dmelFBgn2dmelSymbol[dmelFBgn]
                    f.write(GdxxxID + "\t" + dmelFBgn + "\t" + dmelSymbol + "\t" + "\t".join(expression) + "\n")
                    if not dmelFBgn in dct_ortholog_expression.keys():
                        dct_ortholog_expression[dmelFBgn] = dict()
                    dct_ortholog_expression[dmelFBgn][dxxx] = expression

            """ missed ortholog pairs """
            for dxxxFBgn in exp.oldgeneid2expression.keys():
                expression = exp.oldgeneid2expression[dxxxFBgn]
                if dxxxFBgn in ann.lost1_dxxxFBgn:
                    dmelFBgn = ann.dxxxFBgn2dmelFBgn[dxxxFBgn]
                    dmelSymbol = dmelFBgn2dmelSymbol[dmelFBgn]
                    f.write(dxxxFBgn + "\t" + dmelFBgn + "\t" + dmelSymbol + "\t" + "\t".join(expression) + "\n")
                    if not dmelFBgn in dct_ortholog_expression.keys():
                        dct_ortholog_expression[dmelFBgn] = dict()
                    dct_ortholog_expression[dmelFBgn][dxxx] = expression

    """ merged matrix with the same set of ortholog id """
    with open("../data/expression/merged.121ortholog.genic.nrc.txt", 'w') as f:
        """ header """
        f.write("dmelFBgn" + "\t" + "dmelSymbol")
        for dxxx in sharedinfo.ordered_species:
            f.write("\t" + dct_ortholog_header[dxxx].rstrip())
        f.write("\n")

        """ expression """
        for dmelFBgn in dct_ortholog_expression.keys():
            dmelSymbol = dmelFBgn2dmelSymbol[dmelFBgn]
            if len(dct_ortholog_expression[dmelFBgn].keys()) == 9:
                f.write(dmelFBgn + "\t" + dmelSymbol + "\t")
                for dxxx in sharedinfo.ordered_species:
                    f.write("\t".join(dct_ortholog_expression[dmelFBgn][dxxx]))
                    f.write("\t")
                f.write("\n")



if __name__ == '__main__':
    main()

