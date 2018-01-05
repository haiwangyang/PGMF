#!/usr/bin/env python

"""
Purpose:
    Handling splicing event of spanki

"""
import focalannotation
import focalgene
import sharedinfo
from sharedinfo import exist_file, get_lines

class FocalEvent:
    """FocalEvent object"""
    def __init__(self, species, sex, tissue):
        self.species = species
        self.sex = sex
        self.tissue = tissue
        self.filename = self.species + "_" + sex + "_" + tissue + ".merged.events"
        self.lines = get_lines("../data/event", self.filename)
        self.eventinfo = self.get_eventinfo()

    def get_eventinfo(self):
        """ get eventtion info """
        dct = dict()
        for line in self.lines[1:]:
            (event_id, gene_id, gene_name, eventcode, structure, transcript_id, joinstring, confounding_inccov, confounding_exccov, opt_calc_type, inccov, exccov, inc_sites, exc_sites, inc, exc, totalcov, avgcovpersite, psi) = line.rstrip().split("\t")
            if not gene_id in dct.keys():
                dct[gene_id] = dict()
            dct[gene_id][event_id] = [inc, exc, psi]
        return dct

if __name__ == '__main__':
    #tra_dmelFBgn = sharedinfo.tra_refgeneid
    species = "dyak"
    geneid = sharedinfo.tra_species2geneid[species]
    gen = focalgene.FocalGene(species, geneid, "MDADSSVA")
    eve = FocalEvent(species, "f", "ac")
    eve.eventinfo[geneid]



