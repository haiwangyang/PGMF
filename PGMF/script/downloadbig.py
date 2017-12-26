import urllib.request
import os

def fetch_big_file_from_helix_ftp(filename):
    """ get data of big file from helix ftp
        filename e.g., dyak.fasta, dyak.gtf
    """
    data = urllib.request.urlopen("ftp://helix.nih.gov/pub/haiwang/" + filename)
    if filename.endswith("fasta"):
        folder = "genome"
    elif filename.endswith("gtf"):
        folder = "annotation"
    
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
    print("Trying fetch dyak.fasta from ftp.helix and put it in folder")
    fetch_big_file_from_helix_ftp("dper.SVGpredAdded.v2.gtf")
