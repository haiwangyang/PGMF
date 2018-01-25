#!/usr/bin/env python

"""
Purpose:
    Handling intersectbed (-wo) report of two bed12
    first bed12 is PacBio Iso-seq bed (converted from bam)
    second bed12 is Gtf bed (converted from gtf or v3.gtf)

    gtf versionA is FlyBase gtf
    gtf versionB is updated v3 gtf
"""
from sharedinfo import exist_file, get_lines

def sum_comma_sep_str(string):
    """ return sum value of the comma separated string """
    sum = 0
    for i in string.split(","):
        if not i == "":
            #print(i)
            sum += int(i)
    return sum

def get_id2exonlen(lines):
    """ id can be pacbio isoseqid or gtf transid """
    id2exonlen = dict()
    for line in lines:
        (chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts) = line.rstrip().split("\t")        
        id2exonlen[name] = sum_comma_sep_str(blockSizes)
    return(id2exonlen)

range1 = lambda start, end: range(start, end+1)

def get_exons(chromStart, chromEnd, blockSizes, blockStarts):
    """ parse info from bed12 to get exon blocks """
    blockSizes = [int(i) for i in blockSizes.split(",") if not i == "" ]
    blockStarts = [int(i) for i in blockStarts.split(",") if not i == "" ]
    n = len(blockSizes)
    exons = []
    print("block: " + str(n))
    print(blockSizes,  blockStarts)
    for i in range(n):
        print(i)
        blockStart = blockStarts[i]
        blockSize = blockSizes[i]
        exonStart = chromStart + blockStart
        exonEnd = exonStart + blockSize
        exons.append([exonStart, exonEnd])
    return(exons)

class FocalIntersect:
    """FocalIntersect object"""
    def __init__(self, species, sex, tissue, replicate):
        self.species = species
        self.sex = sex
        self.tissue = tissue
        self.replicate = replicate
        self.name = species + "_" + sex + "_" + tissue + "_" + replicate
        self.filenamePacBioBed = self.name + ".bam.bed"
        print("get lines of " + self.filenamePacBioBed)
        self.linesPacBioBed = get_lines("../data/pacbio", self.filenamePacBioBed)
        print("get isoseqid to exonlen ...")
        self.isoseqid2exonlen = get_id2exonlen(self.linesPacBioBed)        

        self.filenameAnnotationBedA = self.species + ".genePred.bed"
        print("get lines of " + self.filenameAnnotationBedA)
        self.linesAnnotationBedA = get_lines("../data/annotation", self.filenameAnnotationBedA)
        print("get gtf's transid to exonlen ...")
        self.transidA2exonlen = get_id2exonlen(self.linesAnnotationBedA)
        
        self.filenameAnnotationBedB = self.species + ".v3.genePred.bed"
        print("get lines of " + self.filenameAnnotationBedB)
        self.linesAnnotationBedB = get_lines("../data/annotation", self.filenameAnnotationBedB)
        print("get v3.gtf's transid to exonlen ...")
        self.transidB2exonlen = get_id2exonlen(self.linesAnnotationBedB) 
        
        self.filenameA = self.name + ".bam.bed.intersect_gtf"
        self.filenameB = self.name + ".bam.bed.intersect_v3gtf"
        print("get lines of " + self.filenameA)
        self.linesA = get_lines("../data/pacbio", self.filenameA)
        print("get lines of " + self.filenameB)
        self.linesB = get_lines("../data/pacbio", self.filenameB)
        print("get intersection A info ...")
        (self.isoseqid2transidA, self.transidA2isoseqid, self.isoseqidtransidA2info) = self.get_intersectinfo("A")
        print("get intersection B info ...")
        (self.isoseqid2transidB, self.transidB2isoseqid, self.isoseqidtransidB2info) = self.get_intersectinfo("B")

        self.isoseqid2besttransidA = self.get_isoseqid2besttransid("A")
        self.isoseqid2besttransidB = self.get_isoseqid2besttransid("B")

    def get_intersectinfo(self, gtf_version):
        if gtf_version == "A":
            lines = self.linesA
        elif gtf_version == "B":
            lines = self.linesB
        isoseqid2transid = dict()
        transid2isoseqid = dict()
        isoseqid2exonlen = dict()
        transid2exonlen = dict()
        isoseqidtransid2info = dict()

        for line in lines:
            (chrom1, chromStart1, chromEnd1, isoseqid, score1, strand1, thickStart1, thickEnd1, itemRgb1, blockCount1, blockSizes1, blockStarts1, chrom2, chromStart2, chromEnd2, transid, score2, strand2, thickStart2, thickEnd2, itemRgb2, blockCount2, blockSizes2, blockStarts2, intersection)  = line.rstrip().split("\t")
            intersection = int(intersection)
            if not isoseqid in isoseqid2transid.keys():
                isoseqid2transid[isoseqid] = set()
            isoseqid2transid[isoseqid].add(transid)
           
            if not transid in transid2isoseqid.keys():
                transid2isoseqid[transid] = set()
            transid2isoseqid[transid].add(isoseqid)
           
            isoseqidtransid2info[isoseqid + "." + transid] = line
        return(isoseqid2transid, transid2isoseqid, isoseqidtransid2info)

    def get_isoseqid2besttransid(self, gtf_version):
        """ get the transid with the largest overlap """
        if gtf_version == "A":
            isoseqid2transid = self.isoseqid2transidA
            transid2exonlen = self.transidA2exonlen
        elif gtf_version == "B":
            isoseqid2transid = self.isoseqid2transidB
            transid2exonlen = self.transidB2exonlen

        isoseqid2besttransid = dict()
        for isoseqid in isoseqid2transid.keys():
            transids = isoseqid2transid[isoseqid]
            this_transid2exonlen = dict()
            for transid in transids:
                exonlen = int(transid2exonlen[transid])             
                this_transid2exonlen[transid] = exonlen
            inverse = [(value, key) for key, value in this_transid2exonlen.items()]
            isoseqid2besttransid[isoseqid] = max(inverse)[1]         
        return(isoseqid2besttransid)

if __name__ == '__main__':
    ins = FocalIntersect("dmel", "f", "wb", "r1")
    for isoseqid in ins.isoseqid2besttransidA.keys():
        transidA = ins.isoseqid2besttransidA[isoseqid]
        line = ins.isoseqidtransidA2info[isoseqid + "." + transidA]
        # print(isoseqid, transidA, line)
        (chrom1, chromStart1, chromEnd1, isoseqid, score1, strand1, thickStart1, thickEnd1, itemRgb1, blockCount1, blockSizes1, blockStarts1, chrom2, chromStart2, chromEnd2, transid, score2, strand2, thickStart2, thickEnd2, itemRgb2, blockCount2, blockSizes2, blockStarts2, intersection)  = line.rstrip().split("\t")
        chromStart1 = int(chromStart1)
        chromEnd1 = int(chromEnd1)
        chromStart2 = int(chromStart2)
        chromEnd2 = int(chromEnd2)
        exons1 = get_exons(chromStart1, chromEnd1, blockSizes1, blockStarts1)
        #exons2 = get_exons(chromStart2, chromEnd2, blockSizes2, blockStarts2)
        #minStart =  min(chromStart1, chromStart2)
        #maxEnd = max(chromEnd1, chromEnd2)
        print(chromStart1, chromEnd1, blockSizes1, blockStarts1, exons1)
        #for i in range1(minStart, maxEnd):
            
