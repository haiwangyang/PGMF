#!/usr/bin/env python

"""
Purpose:
    Provide supporting information, such as
        functions
        species
        important genes
"""
import downloadbig
from pathlib import Path




""" shared functions """ 

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

def combination_of_two_lists(lst1, lst2):
    """
        ["f", "m"] and ["wb", "go", "re"] => ["f_wb", "f_go", "f_re", "m_wb", "m_go", "m_re"]
    """
    return [e1 + "_" + e2 for e1 in lst1 for e2 in lst2]

""" shared species and gene information """

ordered_species = ['dmel', 'dyak', 'dana', 'dpse', 'dper', 'dwil', 'dmoj', 'dvir', 'dgri']
ordered_tissue8 = ['wb', 'go', 're', 'ge', 'tx', 'dg', 'hd', 'ac']
ordered_tissue7 = ['wb', 'go', 're', 'tx', 'dg', 'hd', 'ac']
ordered_sex = ['f', 'm']

ordered_sexedtissue7 = combination_of_two_lists(ordered_sex, ordered_tissue7)

ordered_sexedtissue8 = combination_of_two_lists(ordered_sex, ordered_tissue8)


tra_refgeneid = "FBgn0003741"
tra_geneid = ["MFBST.6910", "MFBST.7004", "MFBST.9251", "MFBST.13092", "MFBST.6354", "MFBST.7187", "MFBST.10183", "MFBST.10131", "MFBST.4816"]
tra_species2geneid = dict(zip(ordered_species, tra_geneid))

dsx_refgeneid = "FBgn0000504"
dsx_geneid = ['MFBST.8049', 'MFBST.8207', 'MFBST.11795', 'MFBST.1800', 'MFBST.313', 'MFBST.11568', 'MFBST.9018', 'MFBST.2242', 'MFBST.3971' ]
dsx_species2geneid = dict(zip(ordered_species, dsx_geneid))

fru_refgeneid = "FBgn0004652"
fru_geneid = ['FBgn0004652', 'MFBST.7875', 'MFBST.11171', 'MFBST.1237', 'MFBST.847', 'MFBST.11356', 'MFBST.7842', 'MFBST.1901', 'MFBST.2717']
fru_species2geneid = dict(zip(ordered_species, fru_geneid))




""" temperary file locations """

genomefilebiowulf = "/data/yangh13/python/packages/CRoS/CRoS/species_FB2017_03/fa/dxxx.fasta"
annotationfilebiowulf = "/data/yangh13/python/packages/CRoS/CRoS/species_FB2017_03/gtf/dxxx.SVGpredAdded.v2.gtf"
expressionfilebiowulf = "/data/yangh13/python/packages/CRoS/CRoS/htseq/FB2017_03_v2/dxxx.expression.nrc.tab" 

