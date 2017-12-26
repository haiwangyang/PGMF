#!/usr/bin/env python

"""
Purpose:
Obtain features of focal gene

Parameters:
--species (-s): species
--geneid (-g): gene id
--startpeptide (-M): start peptide sequence
--startpeptidelike (-Ml): start peptide sequence like (allow approximation)

Example command:
python3 -i focalgene.py -s dper -g MFBST.6354 -Ml MDADSS

Example command (interactive):
tra = focalgene("dper", "MFBST.6354", "MDADSS")
tra.isoform['FBtr0183594']["exon"][1]
"""

import re
import sys
import os
import argparse
import urllib.request
from pyfaidx import Fasta
from Bio.Seq import Seq

def get_revcom_dna(dna):
    return(str(Seq(dna).reverse_complement()))

def get_pep_and_leftover_from_dna(dna, phase_left):
    """ translate dna into peptide from phase
        phase = 0, from the 1st dna base
        phase = 1, from the 2nd dna base
        phase = 2, from the 3rd dna base
        phased peptide and leftover dna on the left & right are returned
    """
    dna_leftover_left = dna[:phase_left]
    phase_right = len(dna[phase_left:]) % 3
    if phase_right == 0:
        dna_leftover_right = ""
        dna_phased = dna[phase_left:]
    else:
        dna_leftover_right = dna[-phase_right:]
        dna_phased = dna[phase_left:-phase_right]
    pep_phased = str(Seq(dna_phased).translate())
    return(pep_phased, dna_leftover_left, dna_leftover_right)

def get_phase_left_without_stop_codon(dct):
    """ get the phase on the left with correct frame """
    correct_ones = list()
    for phase_left in range(0,3):
        pep_phased, dna_leftover_left, dna_leftover_right = dct[phase_left]
        if not "*" in pep_phased:
            correct_ones.append(phase_left)
    return correct_ones

def get_phase_left_with_startpeptide(dct, startpeptide):
    """ get the phase on the left with correct start peptide"""
    correct_ones = list()
    for phase_left in range(0,3):
        pep_phased, dna_leftover_left, dna_leftover_right = dct[phase_left]
        if startpeptide in pep_phased:
            correct_ones.append(phase_left)
    return correct_ones

def fetch_from_helix_ftp(filename):
    """ get data of big file from helix ftp """
    #data = urllib.request.urlopen("ftp://helix.nih.gov/pub/haiwang/dmel.fasta")
    data = urllib.request.urlopen("ftp://helix.nih.gov/pub/haiwang/" + filename)
    return data


class FocalGene:
    """focalgene object"""
    def __init__(self, species, geneid, startpeptide):
        self.species = species
        self.geneid = geneid
        self.startpeptide = startpeptide
        self.annotation = self.get_annotation_gene()
        self.scaffold = self.get_scaffold()
        self.scaffoldseq = self.get_scaffoldseq()
        self.isoform = self.get_isoform()
        self.isoform = self.update_isoform()
        
    def get_annotation_all(self):
        """ Open annotaiton file """
        with open("../data/annotation/" + self.species + ".SVGpredAdded.v2.gtf", 'r') as f:
            lines = f.readlines()
            return lines

    def get_annotation_gene(self):
        """ Get annotation of focalgene """
        lines_of_focalgene = []
        for line in self.get_annotation_all():
            (scaffold, tag, feature, start, end, scoare, strand, dot, others) = line.rstrip().split("\t")
            this_geneid = re.search('gene_id "(.+?)"', others).group(1)
            if self.geneid == this_geneid:
                lines_of_focalgene.append(line)
        return lines_of_focalgene	

    def get_scaffold(self):
        """ Get scaffold of focalgene """
        return self.annotation[0].split("\t")[0]

    def get_scaffoldseq(self):
        """ Get chromosome sequence of the focal gene """
        return Fasta("../data/genome/" + self.species + ".fasta", as_raw = True)[self.scaffold]    
  
    def get_isoform(self):
        """ Get raw isoform information """
        transid2info = dict()
        for line in self.annotation:
            (scaffold, tag, feature, start, end, scoare, strand, dot, others) = line.rstrip().split("\t")
            this_transid = re.search('transcript_id "(.+?)"', others).group(1)
            if feature == "exon":
                if not this_transid in transid2info.keys():
                    transid2info[this_transid] = dict()
                    transid2info[this_transid]["strand"] = strand
                    transid2info[this_transid]["exonrange"] = set()
                start = int(start)
                end = int(end)
                transid2info[this_transid]["exonrange"].add((start, end))
        return transid2info

    def update_isoform(self):
        """ Update isoform information:
            (1) reorder exon ATG to TAA
            (2) index exon
            (3) fetch sequence
        """
        transid2updatedinfo = dict()
        for transid in self.isoform.keys():
            transid2updatedinfo[transid] = dict()
            transid2updatedinfo[transid]["exon"] = dict()   
            strand = self.isoform[transid]["strand"]
            transid2updatedinfo[transid]["strand"] = strand
            sorted_exonrange = list(self.isoform[transid]["exonrange"])
            if strand == "+":
                sorted_exonrange.sort(key = lambda x: x[0])
            elif strand == "-":
                sorted_exonrange.sort(key = lambda x: x[1], reverse = True)

            exonindex = 1
            isoformdnaseq = ''
            for s, e in sorted_exonrange:
                transid2updatedinfo[transid]["exon"][exonindex] = dict()
                transid2updatedinfo[transid]["exon"][exonindex]["range"] = [s, e]      
                if strand == "+":
                    exonseq = self.scaffoldseq[s-1:e]
                elif strand == "-":
                    exonseq = get_revcom_dna(self.scaffoldseq[s-1:e])
                transid2updatedinfo[transid]["exon"][exonindex]["dnaseq"] = exonseq
                transid2updatedinfo[transid]["exon"][exonindex]["pepseq"] = dict()
                for phase_left in range(0,3):
                    transid2updatedinfo[transid]["exon"][exonindex]["pepseq"][phase_left] = get_pep_and_leftover_from_dna(exonseq, phase_left)
                isoformdnaseq = isoformdnaseq + exonseq
                exonindex += 1
            transid2updatedinfo[transid]["isoformdnaseq"] = isoformdnaseq
        return transid2updatedinfo

if __name__ == '__main__':
    tra = FocalGene("dper", "MFBST.6354", "MDADSSVA")

    # functional tra isoform
    traf = tra.isoform['MFBST.6354.1']
    isoformdanseqf = traf['isoformdnaseq']

    # the first exon
    traf1 = traf['exon'][1]
    phasef1 = get_phase_left_with_startpeptide(traf1["pepseq"], tra.startpeptide)[0]

    print(">dper_tra_functional\n" + get_pep_and_leftover_from_dna(isoformdanseqf, phasef1)[0] + "\n")



    #traf_2 = traf['exon'][2]
    #print(get_phase_left_without_stop_codon(traf_2["pepseq"]))

    #traf_3 = traf['exon'][3]




         
