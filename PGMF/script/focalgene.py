#!/usr/bin/env python

"""
Obtain features of focal gene

Argv:  geneid
"""

import re
import sys
import os
import argparse

parser = argparse.ArgumentParser(description='please provide species and geneid')
parser.add_argument('-s', '--species', type=str)
parser.add_argument('-g', '--geneid', type=str)
args = parser.parse_args()

class focalgene:
    """focalgene object"""
    def __init__(self, species, geneid):
        self.species = species
        self.geneid = geneid

    def get_annotation_all(self):
        """ Open annotaiton file """
        with open("../data/annotation/" + self.species + ".tra.gtf", 'r') as f:
            file_contents = f.read()
            return file_contents

    def get_annotation_gene(self):
        """ Get annotation of gene """
        for line in self.get_annotation_all().split():
            print(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='please provide species and geneid')
    parser.add_argument('-s', '--species', type=str)
    parser.add_argument('-g', '--geneid', type=str)
    args = parser.parse_args()
    tra = focalgene(args.species, args.geneid)
    tra.get_annotation_gene()

