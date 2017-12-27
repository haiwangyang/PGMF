#!/usr/bin/env python

"""
Purpose:
Provide supporting information, such as
species
important genes
"""

ordered_species = ['dmel', 'dyak', 'dana', 'dpse', 'dper', 'dwil', 'dmoj', 'dvir', 'dgri']

dsx_refgeneid = "FBgn0000504"
dsx_geneid = ['MFBST.8049', 'MFBST.8207', 'MFBST.11795', 'MFBST.1800', 'MFBST.313', 'MFBST.11568', 'MFBST.9018', 'MFBST.2242', 'MFBST.3971' ]
dsx_species2geneid = dict(zip(ordered_species, dsx_geneid))

fru_refgeneid = "FBgn0004652"
fru_geneid = ['FBgn0004652', 'MFBST.7875', 'MFBST.11171', 'MFBST.1237', 'MFBST.847', 'MFBST.11356', 'MFBST.7842', 'MFBST.1901', 'MFBST.2717']
fru_species2geneid = dict(zip(ordered_species, fru_geneid))


class SharedInfo:
    """ SharedInfo object """
    def __init__(self):
        self.species = ['dmel', 'dyak', 'dana', 'dpse', 'dper', 'dwil', 'dmoj', 'dvir', 'dgri']
        self.genomefilebiowulf = "/data/yangh13/python/packages/CRoS/CRoS/species_FB2017_03/fa/dxxx.fasta"
        self.annotationfilebiowulf = "/data/yangh13/python/packages/CRoS/CRoS/species_FB2017_03/gtf/dxxx.SVGpredAdded.v2.gtf"
        self.expressionfilebiowulf = "/data/yangh13/python/packages/CRoS/CRoS/htseq/FB2017_03_v2/dxxx.expression.nrc.tab" 

