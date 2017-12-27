#!/usr/bin/env python

"""
Purpose:
    Download big files online
    Files could be dper.fasta or dper.SVGpredAdded.v2.gtf

Example (interactive):
    fetch_big_file_from_helix_ftp("dper.genome")
    fetch_big_file_from_helix_ftp("dper.SVGpredAdded.v2.gtf")
"""

import urllib.request
import os

def fetch_big_file_from_helix_ftp(folder, filename):
    """ get data of big file from helix ftp
        filename e.g., dyak.fasta, dyak.gtf, dyak.expression.nrc.tab
    """
    data = urllib.request.urlopen("ftp://helix.nih.gov/pub/haiwang/" + filename)
    #if filename.endswith("fasta"):
    #    folder = "genome"
    #elif filename.endswith("gtf"):
    #    folder = "annotation"
    #elif filename.endswith("expression.nrc.tab"):
    #    folder = "expression"
    
    # check if the output dir exist, if not create it    
    outputdir = "../data/" + folder
    try:
        os.stat(outputdir)
    except:
        os.mkdir(outputdir)

    with open(outputdir + "/" + filename, 'w', encoding='utf-8') as f:
        for line in data:
            f.write(line.decode('UTF-8'))

if __name__ == '__main__':
    fetch_big_file_from_helix_ftp("dper.fasta")
    fetch_big_file_from_helix_ftp("dper.SVGpredAdded.v2.gtf")
    fetch_big_file_from_helix_ftp("dper.expression.nrc.tab")
