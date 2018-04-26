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

def get_elements(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
        elements = list()
        for line in lines:
             elements.append(line.rstrip())
        return(elements)


def get_A2B(filepath):
    """
       input a text table sep by "\t"
       return a dict with
       col1 key
       col2 value
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
        dct = dict()
        for line in lines:
             elements = line.rstrip().split("\t")
             dct[elements[0]] = elements[1]
        return(dct)

def get_M2N(filepath, m, n):
    """
       input a text table sep by "\t"
       return a dict with
       m: key col index
       n: value col index
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
        dct = dict()
        for line in lines:
             elements = line.rstrip().split("\t")
             dct[elements[m]] = elements[n]
        return(dct)


def combination_of_two_lists(lst1, lst2):
    """
        ["f", "m"] and ["wb", "go", "re"] => ["f_wb", "f_go", "f_re", "m_wb", "m_go", "m_re"]
    """
    return [e1 + "_" + e2 for e1 in lst1 for e2 in lst2]

def DxxxGID_to_YOgnID(DxxxGID):
    """ convert DxxxGID to YOgnID """
    return("YOgn" + DxxxGID[1:3].upper() + DxxxGID[-5:])

def DxxxTID_to_YOtrID(DxxxTID):
    """ convert DxxxTID to YOtrID """
    return("YOtr" + DxxxTID[1:3].upper() + DxxxTID[-6:])

def jaccard(exonmap1, exonmap2):
    """ get jaccard between exon map
        1.3 and 2.4
        jaccard = 2/4 = 0.5
    """
    union_sum = 0
    intersection_sum = 0

    dct1 = dict()
    for se in exonmap1:
        s, e = se.split(".")
        for i in range(int(s), int(e) + 1):
            if not i in dct1.keys():
                dct1[i] = 0
            dct1[i] += 1

    dct2 = dict()
    for se in exonmap2:
        s, e = se.split(".")
        for i in range(int(s), int(e) + 1):
            if not i in dct2.keys():
                dct2[i] = 0
            dct2[i] += 1

    st = set()
    for ii in [dct1.keys(), dct2.keys()]:
        for i in ii:
            st.add(i)

    union_sum = len(st)
    for i in st:
        if i in dct1.keys() and i in dct2.keys():
            intersection_sum += 1

    j = intersection_sum / union_sum
    return(j)

def get_GSM():
    """ get sample to GSM and sampleInGEO """
    dct = dict()
    for line in get_lines("../data/GSM", "GSE.txt"):
        GSM, sampleInGEO, sample = line.rstrip().split("\t")
        if not sample in dct.keys():
            dct[sample] = []
        dct[sample].append([GSM, sampleInGEO])
    return(dct)

""" shared species and gene information """
typical_strains = ['w1118', 'dyak']
# typical_strains = ['w1118', 'dyak', 'dana', 'dpse', 'dper', 'dwil', 'dmoj', 'dvir', 'dgriG1']
ordered_species = ['dmel', 'dyak', 'dana', 'dpse', 'dper', 'dwil', 'dmoj', 'dvir', 'dgri']
species2dxxx = get_A2B("../data/species/species2dxxx.txt")

samples = get_elements("../data/sample/sample.list")
sample2GSM = get_GSM()

ordered_tissue8 = ['wb', 'go', 're', 'ge', 'tx', 'dg', 'hd', 'ac']
ordered_tissue7 = ['wb', 'go', 're', 'tx', 'dg', 'hd', 'ac']
ordered_sex = ['f', 'm']

ordered_sexedtissue7 = combination_of_two_lists(ordered_sex, ordered_tissue7)

ordered_sexedtissue8 = combination_of_two_lists(ordered_sex, ordered_tissue8)

# dmelFBgnID2dmelSymbol = FocalAnnotation("dmel").dmelFBgn2dmelSymbol


tra_refgeneid = "FBgn0003741"
tra_geneid = ["MFBST.6910", "MFBST.7004", "MFBST.9251", "MFBST.13092", "MFBST.6354", "MFBST.7187", "MFBST.10183", "MFBST.10131", "MFBST.4816"]
tra_species2geneid = dict(zip(ordered_species, tra_geneid))

dsx_refgeneid = "FBgn0000504"
dsx_geneid = ['MFBST.8049', 'MFBST.8207', 'MFBST.11795', 'MFBST.1800', 'MFBST.313', 'MFBST.11568', 'MFBST.9018', 'MFBST.2242', 'MFBST.3971' ]
dsx_species2geneid = dict(zip(ordered_species, dsx_geneid))

fru_refgeneid = "FBgn0004652"
fru_geneid = ['FBgn0004652', 'MFBST.7875', 'MFBST.11171', 'MFBST.1237', 'MFBST.847', 'MFBST.11356', 'MFBST.7842', 'MFBST.1901', 'MFBST.2717']
fru_species2geneid = dict(zip(ordered_species, fru_geneid))

pacbio_sample = ["dgri_f_wb_r1", "dgri_f_wb_r2", "dgri_m_wb_r1", "dgri_m_wb_r2", "dmel_f_go_r1", "dmel_f_wb_r1", "dmel_m_go_r1", "dmel_m_wb_r1"]


""" temperary file locations """

genomefilebiowulf = "/data/yangh13/python/packages/CRoS/CRoS/species_FB2017_03/fa/dxxx.fasta"
annotationfilebiowulf = "/data/yangh13/python/packages/CRoS/CRoS/species_FB2017_03/gtf/dxxx.SVGpredAdded.v2.gtf"
expressionfilebiowulf = "/data/yangh13/python/packages/CRoS/CRoS/htseq/FB2017_03_v2/dxxx.expression.nrc.tab" 

