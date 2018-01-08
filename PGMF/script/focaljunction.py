#!/usr/bin/env python

"""
Purpose:
    Handling normalized read count of genes

"""
import focalannotation
import focalgene
import sharedinfo
import pandas as pd
from sharedinfo import exist_file, get_lines

def get_junction_of_species_by_partiallocation(species, partiallocation):
    """
       species (e.g., dyak)
       location is partial location (e.g., 3L_17475)
       output is all possible junctions with coverage
    """
    geneid = sharedinfo.tra_species2geneid[species]
    gen = focalgene.FocalGene(species, geneid, "M")
    for sex in sharedinfo.ordered_sex:
        for tissue in sharedinfo.ordered_tissue7:
            jun = FocalJunction(species, sex, tissue)
            for range in jun.juncinfo.keys():
                if range.startswith(partiallocation):
                    print(species + "_" + sex + "_" + tissue, partiallocation, range, jun.juncinfo[range])

def get_junction_of_species_by_location(species, location):
    """
       species (e.g., dyak)
       location is whole range (e.g., 3L_17475357_17475426_+)
       output is a dict of junctions with coverage
    """
    geneid = sharedinfo.tra_species2geneid[species]
    gen = focalgene.FocalGene(species, geneid, "M")
    dct = dict()
    for sex in sharedinfo.ordered_sex:
        for tissue in sharedinfo.ordered_tissue7:
            dct[species + "_" + sex + "_" + tissue] = dict()
            jun = FocalJunction(species, sex, tissue)
            
            if location in jun.juncinfo:
                dct[species + "_" + sex + "_" + tissue][location] = jun.juncinfo[location]
            else:
                dct[species + "_" + sex + "_" + tissue][location] = 0
    return dct

def merge_dcts(dcts):
    mdct = dict()
    for d in dcts:
        for sst in d.keys(): # species_sex_tissue
            if not sst in mdct.keys():
                mdct[sst] = dict()
            for range in d[sst].keys():
                mdct[sst][range] = d[sst][range]
    return mdct

class FocalJunction:
    """FocalJunction object"""
    def __init__(self, species, sex, tissue):
        self.species = species
        self.sex = sex
        self.tissue = tissue
        # self.filename = self.species + "_" + sex + "_" + tissue + ".merged.juncs" # spanki juncs
        self.filename = self.species + "_" + sex + "_" + tissue + ".sorted.junc.bed"
        self.lines = get_lines("../data/junction", self.filename)
        self.juncinfo = self.get_juncinfo()

    def get_juncinfo(self):
        """  get junction info """
        dct = dict()
        for line in self.lines:
              # (juncid, dinucleotide, intron_size, annostatus, gmcode, regcode, geneassign, cov, lirt, rirt, irt, dncov, ancov, numsamps) = line.rstrip().split("\t") # spanki junc
            (chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts) = line.rstrip().split("\t")
            chromStart = int(chromStart)
            chromEnd = int(chromEnd)
            bs1, bs2 = blockSizes.split(",")
            bs1 = int(bs1)
            bs2 = int(bs2)
            juncstart = chromStart + bs1 + 1
            juncend = chromEnd - bs2
            dct[chrom + ":" + str(juncstart) + "-" + str(juncend) + "_" + strand] = score
        return dct

if __name__ == '__main__':
    
    # consensus splicing junction of tra in dyak
    dct_cs = get_junction_of_species_by_location("dyak", "3L:17475357-17475426_+")

    # alternative splicing junction (short) of tra in dyak
    dct_as1 = get_junction_of_species_by_location("dyak", "3L:17474772-17474844_+")

    # alternative splicing junction (long) of tra in dyak
    dct_as2 = get_junction_of_species_by_location("dyak", "3L:17474772-17475015_+")

    mdct = merge_dcts([dct_cs, dct_as1, dct_as2])

    mpd = pd.DataFrame.from_dict(mdct)

    mpd.to_csv("../data/output/dyak.tra.junc.summary.txt", sep="\t")
