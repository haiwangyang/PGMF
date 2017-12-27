#!/usr/bin/env python

"""
Purpose:
Annotation ID

"""
import re
import downloadbig
import sharedinfo
from pathlib import Path

def exist_file(folder, filename):
    """ if the file exist """
    path = Path(folder + "/" + filename)
    try:
        path.resolve()    
    except:
        return False
    else:
        return True

def get_lines(folder, filename):
    """ get lines of the file """
    if not exist_file(folder, filename):
        downloadbig.fetch_big_file_from_helix_ftp(folder, filename)

    with open(folder + "/" + filename, 'r') as f:
        lines = f.readlines()
        return lines

def get_id(others, tag):
    """ get id from gtf others field
        id could be geneid, refgeneid, ...
    """
    try:
        id = re.search(tag + ' "(.+?)"', others).group(1)
    except AttributeError:
        id = ""
    return id
        

class FocalAnnotation:
    """focalgene object"""
    def __init__(self, species):
        self.species = species
        self.filename = species + ".SVGpredAdded.v2.gtf"
        self.lines = get_lines("../data/annotation", self.filename)
        self.geneid, self.refgeneid, self.geneid2refgeneid, self.refgeneid2geneid  = self.get_geneidtable()
        self.geneid2refgeneid121, self.refgeneid2geneid121 = self.get_121_id()

    def get_geneidtable(self):
        """ get geneid, refgeneid, and relationship """
        st_geneid = set()
        st_refgeneid = set()
        dct_geneid2refgeneid = dict()
        dct_refgeneid2geneid = dict()
        for line in self.lines:
            (scaffold, tag, feature, start, end, scoare, strand, dot, others) = line.rstrip().split("\t")
            this_geneid = get_id(others, 'gene_id')
            this_refgeneid = get_id(others, 'ref_gene_id')

            count = 0
            if not this_geneid == "":
                st_geneid.add(this_geneid)
                count += 1
            if not this_refgeneid == "":
                st_refgeneid.add(this_refgeneid)
                count += 1
            if count == 2:
                if not this_geneid in dct_geneid2refgeneid.keys():
                    dct_geneid2refgeneid[this_geneid] = set()
                if not this_refgeneid in dct_refgeneid2geneid.keys():
                    dct_refgeneid2geneid[this_refgeneid] = set()
                
                dct_geneid2refgeneid[this_geneid].add(this_refgeneid)
                dct_refgeneid2geneid[this_refgeneid].add(this_geneid)

        return st_geneid, st_refgeneid, dct_geneid2refgeneid, dct_refgeneid2geneid

    def get_121_id(self):
        """ get one to one (geneid to refgeneid)"""

        geneid_one = set()
        for geneid in self.geneid2refgeneid.keys():
            if len(self.geneid2refgeneid[geneid]) == 1:
                geneid_one.add(geneid)

        refgeneid_one = set()
        for refgeneid in self.refgeneid2geneid.keys():
            if len(self.refgeneid2geneid[refgeneid]) == 1:
                refgeneid_one.add(refgeneid)

        dct_geneid2refgeneid121 = dict()
        dct_refgeneid2geneid121 = dict()
        for geneid in geneid_one:
            for refgeneid in self.geneid2refgeneid[geneid]:
                if refgeneid in refgeneid_one:
                    dct_geneid2refgeneid121[geneid] = refgeneid
                    dct_refgeneid2geneid121[refgeneid] = geneid

        return dct_geneid2refgeneid121, dct_refgeneid2geneid121

if __name__ == '__main__':
#    for dxxx in sharedinfo.ordered_species:
    for dxxx in ["dyak", ]:
        print(dxxx)
        ann = FocalAnnotation(dxxx)


