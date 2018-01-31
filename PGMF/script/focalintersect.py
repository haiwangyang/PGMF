#!/usr/bin/env python

"""
Purpose:
    Handling intersectbed (-wo) report of two bed12
    first bed12 is PacBio Iso-seq bed (converted from bam)
    second bed12 is Gtf bed (converted from gtf or v3.gtf)
    
    gtf versionA is FlyBase gtf
    gtf versionB is updated v3 gtf

    potential interesting region
    dmel X:11165258-11167748
"""

from __future__ import print_function, division
from sharedinfo import exist_file, get_lines
import sharedinfo
import sys

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

""" new range1 """
range1 = lambda start, end: range(start, end+1)


def inside_exons(i, exons):
    """ check if position i inside structure of exons """
    for exon in exons:
        if i in range1(exon[0], exon[1]):
            return(True)
    return(False)

def get_exons(chromStart, chromEnd, blockSizes, blockStarts):
    """ parse info from bed12 to get exon blocks """
    blockSizes = [int(i) for i in blockSizes.split(",") if not i == "" ]
    blockStarts = [int(i) for i in blockStarts.split(",") if not i == "" ]
    n = len(blockSizes)
    exons = []
    #print("block: " + str(n))
    #print(blockSizes,  blockStarts)
    for i in range(n):
        #print(i)
        blockStart = blockStarts[i]
        blockSize = blockSizes[i]
        exonStart = chromStart + blockStart
        exonEnd = exonStart + blockSize
        exons.append([exonStart, exonEnd])
    return(exons)

""" convert [[1, 2], [3, 4]] to {'3.4', '1.2'} """
convert_2dlst_to_set = lambda lst: set([float(str(i) + "." + str(ii)) for i,ii in lst])
convert_set_to_2dlst = lambda sat: [list(map(int, i.split("."))) for i in sat]

def return_first_and_last_exon(exons):
    """ float format of exons
        {7.8, 5.6, 3.4, 1.2}
        only consider genome coordinate
        return (1.2, 7.8)
    """
    lst_exons = list(exons)
    lst_exons.sort(key=int)
    return(lst_exons[0], lst_exons[-1])    


def report_exon_overlap(strand1, exons1, strand2, exons2):
    """ report summary of exons1 exons2 overlap """
    #print(strand1)
    #print(exons1)
    #print(exons2)
    exons1 = convert_2dlst_to_set(exons1)
    first_exon1, last_exon1 = return_first_and_last_exon(exons1)
    exons2 = convert_2dlst_to_set(exons2)
    first_exon2, last_exon2 = return_first_and_last_exon(exons2)
    
    dct_report = dict()
    if not first_exon1 == first_exon2:
        """ first exon of isoseq and annotated-gene-model are not exactly the same """
        if str(first_exon1).split(".")[1] == str(first_exon2).split(".")[1]:
            """ if first intron-start boundary is the same """
            if int(str(first_exon1).split(".")[0]) > int(str(first_exon2).split(".")[0]):
                """ if isoseq first exon is shorter """
                if strand1 == "+":
                    dct_report[5] = "partial_inside"
                else:
                    dct_report[3] = "partial_inside"
            else:
                """ if isoseq first exon is longer """
                if strand1 == "+":
                    dct_report[5] = "partial_outside"
                else:
                    dct_report[3] = "partial_outside"
        else:
            if strand1 == "+":
                dct_report[5] = "different"
            else:
                dct_report[3] = "different"
    else:
        if strand1 == "+":
            dct_report[5] = "same"
        else:
            dct_report[3] = "same"

    if not last_exon1 == last_exon2:
        """ last exon of isoseq and annotated-gene-model are not exactly the same """
        if str(last_exon1).split(".")[0] == str(last_exon2).split(".")[0]:
            """ if last intron-end boundary is the same """
            if int(str(last_exon1).split(".")[1]) < int(str(last_exon2).split(".")[1]):
                """ if isoseq first exon is shorter """
                if strand1 == "+":
                    dct_report[3] = "partial_inside"
                else:
                    dct_report[5] = "partial_inside"
            else:
                """ if isoseq first exon is longer """
                if strand1 == "+":
                    dct_report[3] = "partial_outside"
                else:
                    dct_report[5] = "partial_outside"
        else:
            if strand1 == "+":
                dct_report[3] = "different"
            else:
                dct_report[5] = "different" 
    else:
        if strand1 == "+":
            dct_report[3] = "same"
        else:
            dct_report[5] = "same"
    return(dct_report[5], dct_report[3])

def summarize_jaccard(sample):
    """ generate summary of jaccard for each isoseq """
    (species, sex, tissue, replicate) = sample.split("_")
    ins = FocalIntersect(species, sex, tissue, replicate)
    with open("../data/output/" + ins.name + ".A.txt", 'w') as f:
        for isoseqid in ins.jaccard_infoA.keys():
             intersection_sum, union_sum, jaccard = ins.jaccard_infoA[isoseqid]
             transid = ins.isoseqid2besttransidA[isoseqid]
             f.write(isoseqid + "\t" + transid + "\t" + str(intersection_sum) + "\t" + str(union_sum) + "\t" + str(jaccard) + "\n")

    with open("../data/output/" + ins.name + ".B.txt", 'w') as f:
        for isoseqid in ins.jaccard_infoB.keys():
             intersection_sum, union_sum, jaccard = ins.jaccard_infoB[isoseqid]
             transid = ins.isoseqid2besttransidB[isoseqid]
             f.write(isoseqid + "\t" + transid + "\t" + str(intersection_sum) + "\t" + str(union_sum) + "\t" + str(jaccard) + "\n")


def summarize_end_for_all_eight_samples():
    """ summarize end isoseq and gene model report """
    for sample in sharedinfo.pacbio_sample:
        (species, sex, tissue, replicate) = sample.split("_")
        ins = FocalIntersect(species, sex, tissue, replicate)
        with open("../data/output/" + ins.name + ".end.A.txt", 'w') as f:
            for isoseqid in ins.end_infoA.keys():
                 exonnum1, exonnum2, five_desc, three_desc = ins.end_infoA[isoseqid]
                 transid = ins.isoseqid2besttransidA[isoseqid]
                 f.write(isoseqid + "\t" + transid + "\t" + exonnum1 + "\t" + exonnum2 + "\t" + five_desc + "\t" + three_desc + "\n")

        with open("../data/output/" + ins.name + ".end.B.txt", 'w') as f:
            for isoseqid in ins.end_infoB.keys():
                 exonnum1, exonnum2, five_desc, three_desc = ins.end_infoB[isoseqid]
                 transid = ins.isoseqid2besttransidB[isoseqid]
                 f.write(isoseqid + "\t" + transid + "\t" + exonnum1 + "\t" + exonnum2 + "\t" + five_desc + "\t" + three_desc + "\n")

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
        
        """ get jaccard info such as intersection union jaccard for each isoseqid """
        print("get jaccard A info ...")
        # self.jaccard_infoA = self.get_jaccard_info("A")
        print("get jaccard B info ...")
        # self.jaccard_infoB = self.get_jaccard_info("B")
        
        """ get end info for each isoseqid """
        print("get end A info ...")
        self.end_infoA = self.get_end_info("A")
        print("get end B info ...")
        self.end_infoB = self.get_end_info("B")


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
                #print(gtf_version, transid, str(exonlen))    
                this_transid2exonlen[transid] = exonlen
            inverse = [(value, key) for key, value in this_transid2exonlen.items()]
            isoseqid2besttransid[isoseqid] = max(inverse)[1]         
        return(isoseqid2besttransid)

    def get_end_info(self, gtf_version):
        """ get the end info for isoseq cmp annotation gene model """
        if gtf_version == "A":
            isoseqid2besttransid = self.isoseqid2besttransidA
            isoseqidtransid2info = self.isoseqidtransidA2info
        elif gtf_version == "B":
            isoseqid2besttransid = self.isoseqid2besttransidB
            isoseqidtransid2info = self.isoseqidtransidB2info

        isoseqid2endinfo = dict()
        count = 0
        for isoseqid in isoseqid2besttransid.keys():
            if count % 1000 == 0:
                print(count)
            count += 1
            transid = isoseqid2besttransid[isoseqid]
            line = isoseqidtransid2info[isoseqid + "." + transid]
            #print("#################################")
            #print(isoseqid, transid, line)
            (chrom1, chromStart1, chromEnd1, isoseqid, score1, strand1, thickStart1, thickEnd1, itemRgb1, blockCount1, blockSizes1, blockStarts1, chrom2, chromStart2, chromEnd2, transid, score2, strand2, thickStart2, thickEnd2, itemRgb2, blockCount2, blockSizes2, blockStarts2, intersection)  = line.rstrip().split("\t")
            chromStart1 = int(chromStart1)
            chromEnd1 = int(chromEnd1)
            chromStart2 = int(chromStart2)
            chromEnd2 = int(chromEnd2)
            exons1 = get_exons(chromStart1, chromEnd1, blockSizes1, blockStarts1)
            exons2 = get_exons(chromStart2, chromEnd2, blockSizes2, blockStarts2)
            minStart =  min(chromStart1, chromStart2)
            maxEnd = max(chromEnd1, chromEnd2)
            five_prime_end, three_prime_end = (report_exon_overlap(strand1, exons1, strand2, exons2))
            isoseqid2endinfo[isoseqid] = [str(len(exons1)), str(len(exons2)), five_prime_end, three_prime_end]
        return(isoseqid2endinfo)

    def get_jaccard_info(self, gtf_version):
        """ get the transid with the largest overlap """
        if gtf_version == "A":
            isoseqid2besttransid = self.isoseqid2besttransidA
            isoseqidtransid2info = self.isoseqidtransidA2info
        elif gtf_version == "B":
            isoseqid2besttransid = self.isoseqid2besttransidB
            isoseqidtransid2info = self.isoseqidtransidB2info

        isoseqid2jaccardinfo = dict()
        count = 0
        for isoseqid in isoseqid2besttransid.keys():
            if count % 1000 == 0:
                print(count)
            count += 1
            transid = isoseqid2besttransid[isoseqid]
            line = isoseqidtransid2info[isoseqid + "." + transid]
            print("#################################")
            print(isoseqid, transid, line)
            (chrom1, chromStart1, chromEnd1, isoseqid, score1, strand1, thickStart1, thickEnd1, itemRgb1, blockCount1, blockSizes1, blockStarts1, chrom2, chromStart2, chromEnd2, transid, score2, strand2, thickStart2, thickEnd2, itemRgb2, blockCount2, blockSizes2, blockStarts2, intersection)  = line.rstrip().split("\t")
            chromStart1 = int(chromStart1)
            chromEnd1 = int(chromEnd1)
            chromStart2 = int(chromStart2)
            chromEnd2 = int(chromEnd2)
            exons1 = get_exons(chromStart1, chromEnd1, blockSizes1, blockStarts1)
            exons2 = get_exons(chromStart2, chromEnd2, blockSizes2, blockStarts2)
            minStart =  min(chromStart1, chromStart2)
            maxEnd = max(chromEnd1, chromEnd2)

            union_sum = 0
            intersection_sum = 0
            for i in range1(minStart, maxEnd):
                first = inside_exons(i, exons1)
                second = inside_exons(i, exons2)
                if first + second > 0:
                    union_sum += 1
                    if first + second == 2:
                        intersection_sum += 1
            jaccard = intersection_sum / union_sum
            isoseqid2jaccardinfo[isoseqid] = [intersection_sum, union_sum, jaccard]
        return(isoseqid2jaccardinfo)

    def generate_jaccard0_isoseq_bed(self):
        """ get the isoseq reads with no overlap with existing gene model in v3 annotation """
        all = set(self.isoseqid2exonlen.keys())
        notwant = set(self.isoseqid2besttransidB.keys())
        want = all - notwant
        want_lines = []
        with open("../data/pacbio/" + self.name + ".B.j0.bed", 'w') as f:
            for line in self.linesPacBioBed:
                (chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts) = line.rstrip().split("\t")
                if name in want:
                    f.write(line)

                


if __name__ == '__main__':
    sample = sys.argv[1] # dmel_f_wb_r1
    # summarize_jaccard(sample)
    summarize_end_for_all_eight_samples()

    (species, sex, tissue, replicate) = sample.split("_")
    ins = FocalIntersect(species, sex, tissue, replicate)
    ins.generate_jaccard0_isoseq_bed()
