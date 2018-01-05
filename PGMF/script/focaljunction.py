#!/usr/bin/env python

"""
Purpose:
    Handling normalized read count of genes

"""
import focalannotation
import focalgene
import sharedinfo
from sharedinfo import exist_file, get_lines

class FocalJunction:
    """FocalJunction object"""
    def __init__(self, species, sex, tissue):
        self.species = species
        self.sex = sex
        self.tissue = tissue
        self.filename = self.species + "_" + sex + "_" + tissue + ".merged.juncs"
        self.lines = get_lines("../data/junction", self.filename)
        self.juncinfo = self.get_juncinfo()

    def get_juncinfo(self):
        """ get junction info """
        dct = dict()
        for line in self.lines:
            (juncid, dinucleotide, intron_size, annostatus, gmcode, regcode, geneassign, cov, lirt, rirt, irt, dncov, ancov, numsamps) = line.rstrip().split("\t")
            if not geneassign in dct.keys():
                dct[geneassign] = dict()
            dct[geneassign][juncid] = cov
        return dct

if __name__ == '__main__':
    species = "dyak"
    geneid = sharedinfo.tra_species2geneid[species]
    gen = focalgene.FocalGene(species, geneid, "MDADSSVA")
    jun = FocalJunction(species, "m", "wb")
    jun.juncinfo[geneid]

