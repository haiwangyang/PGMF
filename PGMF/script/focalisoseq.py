#!/usr/bin/env python

"""
Purpose:
    Handling pacbio isoseq

"""
import sharedinfo
import re
import focalintersect
from pyfaidx import Fasta
import pysam
from collections import Counter
from sharedinfo import exist_file, get_lines

def summarize_polyA(fasta):
    """ summarize polyA type AAAAAAAAA or TTTTTTTTTT or others """
    lst = []
    for name in fasta.keys():
        seq = str(fasta[name])
        L = seq[0:10]
        R = seq[-10:]
        end = L + R
        most_common_char = Counter(end).most_common(1)[0][0]
        Ln = Counter(L)[most_common_char]
        Rn = Counter(R)[most_common_char]
        if Ln > Rn:
            m = re.search('^(' + most_common_char + '+)', seq)
            if m:
                lst.append(["L", most_common_char, m.group(1), name, seq])
            else:
                lst.append(["L", most_common_char, "-", name, seq])
        else:
            m = re.search('(' + most_common_char + '+)$', seq)
            if m:
                lst.append(["R", most_common_char, m.group(1), name, seq])
            else:
                lst.append(["R", most_common_char, "-", name, seq])
    return(lst) 

def printout_polyA_summary(sample):
    """ analyze the list of polyA
    """
    species, sex, tissue, replicate = sample.split("_")
    fi = FocalIsoseq(species, sex, tissue, replicate)
    for i in fi.polyA:
        if i[1] == 'A' or i[1] == 'T':
            print("pass", sample, i[0], i[1], len(i[2]))
        else:
            print("fail", sample, i[0], i[1], len(i[2]))

def seq2polyA(seq):
    """ input seq
        output polyA report
        (1) pass or not
        (2) Left or Right
        (3) most common character
        (4) length of polyA
        (5) length of isoseq
    """
    lst = []
    L = seq[0:10]
    R = seq[-10:]
    end = L + R
    most_common_char = Counter(end).most_common(1)[0][0]
    Ln = Counter(L)[most_common_char]
    Rn = Counter(R)[most_common_char]
    if Ln > Rn:
        m = re.search('^(' + most_common_char + '+)', seq)
        if m:
            lst = ["L", most_common_char, m.group(1)]
        else:
            lst = ["L", most_common_char, "-"]
    else:
        m = re.search('(' + most_common_char + '+)$', seq)
        if m:
            lst = ["R", most_common_char, m.group(1)]
        else:
            lst = ["R", most_common_char, "-"]

    lst2 = []
    if lst[1] == 'A' or lst[1] == 'T':
        lst2 = ["pass", lst[0], lst[1], str(len(lst[2])), str(len(seq))]
    else:
        lst2 = ["fail", lst[0], lst[1], str(len(lst[2])), str(len(seq))]
    return lst2

def printout_polyA_len_for_sample(sample):
    """ print out polyA len and isoseq len for different jaccard type """
    species, sex, tissue, replicate = sample.split("_")
    fi = FocalIsoseq(species, sex, tissue, replicate)
    for isoseqid in fi.jaccard_zero_isoseqid:
        seq = str(fi.fasta[isoseqid])
        intronnum = fi.isoseqid2intronnum[isoseqid]
        print("jaccard==0" + "\t" + str(intronnum) + "\t" + "\t".join(seq2polyA(seq)))

    for isoseqid in fi.jaccard_plus_isoseqid:
        seq = str(fi.fasta[isoseqid])
        intronnum = fi.isoseqid2intronnum[isoseqid]
        print("jaccard>0" + "\t" + str(intronnum) + "\t" + "\t".join(seq2polyA(seq)))

    for isoseqid in fi.unmapped_isoseqid:
        seq = str(fi.fasta[isoseqid])
        print("unmapped" + "\t" + "NA" + "\t" + "\t".join(seq2polyA(seq)))

def generate_jaccard_type_isoseqid_bed_for_sample(sample):
    """ print out isoseq bed for different jaccard type """
    species, sex, tissue, replicate = sample.split("_")
    fi = FocalIsoseq(species, sex, tissue, replicate)
    with open("../data/pacbio/jaccard_zero." + sample + ".bed", "w") as f:
        for isoseqid in fi.jaccard_zero_isoseqid:
            f.write(fi.isoseqid2bambedline[isoseqid])

    with open("../data/pacbio/jaccard_plus." + sample + ".bed", "w") as f:
        for isoseqid in fi.jaccard_plus_isoseqid:
            f.write(fi.isoseqid2bambedline[isoseqid])


def print_unique_fasta_number():
    """ calculate unique fasta, and remove redundancy """
    st = set()
    dgri_samples = ["dgri_f_wb_r1", "dgri_m_wb_r1", "dgri_f_wb_r2", "dgri_m_wb_r2"]
    dmel_samples =  ["dmel_f_go_r1", "dmel_m_go_r1", "dmel_f_wb_r1", "dmel_m_wb_r1"]
    samples = {"dgri":dgri_samples, "dmel":dmel_samples}
    for dxxx in ['dmel', 'dgri']:
        st = set()
        for sample in samples[dxxx]:
            fasta = Fasta("../data/pacbio/" +  sample + ".fasta")
            for name in fasta.keys():
                seq = fasta[name]
                st.add(seq)
        print(len(st))



class FocalIsoseq:
    """FocalIsoseq object"""
    def __init__(self, species, sex, tissue, replicate):
        self.species = species
        self.sex = sex
        self.tissue = tissue
        self.replicate = replicate
        self.sample = species + "_" + sex + "_" + tissue + "_" + replicate
        self.fasta = Fasta("../data/pacbio/" +  self.sample + ".fasta")
        self.bam = pysam.AlignmentFile("../data/pacbio/" + self.sample + ".bam", "rb")
        self.polyA = summarize_polyA(self.fasta)
        self.all_isoseqid = set(self.fasta.keys())
        self.mapped_isoseqid = self.get_mapped_isoseqid()
        self.unmapped_isoseqid = self.all_isoseqid - self.mapped_isoseqid
        self.jaccard_plus_isoseqid = self.get_jaccard_plus_isoseqid()
        self.jaccard_zero_isoseqid = self.mapped_isoseqid - self.jaccard_plus_isoseqid
        self.isoseqid2bambedline = self.get_isoseqid2bambedline()
        self.isoseqid2intronnum = self.get_isoseqid2intronnum()
        self.isoseqid2position = self.get_isoseqid2position()

    def get_jaccard_plus_isoseqid(self):
        """ get jaccard > 0 isoseqid """
        lines = get_lines("../data/output", self.species + "_" + self.sex + "_" + self.tissue + "_" + self.replicate + ".B.txt")
        st = set()
        for line in lines:
            (position_isoseqid, transid, intersection, union, jaccard) = line.rstrip().split("\t")
            position, isoseqid = position_isoseqid.split(".")
            if float(jaccard) > 0:
                st.add(isoseqid)
        return st

    def get_mapped_isoseqid(self):
        """ 
            return isoseqid in bam file
        """
        st = set()
        for read in self.bam:
            isoseqid = read.qname
            st.add(isoseqid)
        return st

    def get_isoseqid2bambedline(self):
        lines = get_lines("../data/pacbio", self.sample + ".bam.bed")
        dct = dict()
        for line in lines:
            (chrom, chromStart, chromEnd, position_isoseqid, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts) = line.rstrip().split("\t")
            position, isoseqid = position_isoseqid.split(".")
            dct[isoseqid] = line
        return dct            

    def get_isoseqid2position(self):
        lines = get_lines("../data/pacbio", self.sample + ".bam.bed")
        dct = dict()
        for line in lines:
            (chrom, chromStart, chromEnd, position_isoseqid, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts) = line.rstrip().split("\t")
            position, isoseqid = position_isoseqid.split(".")
            dct[isoseqid] = position
        return dct

    def get_isoseqid2intronnum(self):
        """ get intron information from bam.bed """
        lines = get_lines("../data/pacbio", self.sample + ".bam.bed")
        dct = dict()
        for line in lines:
            (chrom, chromStart, chromEnd, position_isoseqid, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts) = line.rstrip().split("\t") 
            position, isoseqid = position_isoseqid.split(".")
            exons = focalintersect.get_exons(int(chromStart), int(chromEnd), blockSizes, blockStarts)
            num_exons = len(exons)
            dct[isoseqid] = num_exons - 1
        return dct

def main():
    #""" printout polyA len """
    #for sample in sharedinfo.pacbio_sample:
    #    printout_polyA_summary(sample)

    """ check len distribution of polyA and isoseq """
    # printout_polyA_len_for_sample("dmel_m_go_r1")

    for sample in sharedinfo.pacbio_sample:
        generate_jaccard_type_isoseqid_bed_for_sample(sample)

if __name__ == '__main__':
    # main()
    # fi = FocalIsoseq("dgri", "f", "wb", "r1")
    print_unique_fasta_number()


