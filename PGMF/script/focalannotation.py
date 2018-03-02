#!/usr/bin/env python

"""
Purpose:
    Handling annotation ID
    ortholog one2one (old and new) also included in the FocalAnnotation object
"""
import re
import downloadbig
import sharedinfo
from sharedinfo import exist_file, get_lines

def get_id(others, tag):
    """ get id from gtf others field
        id could be geneid, refgeneid, ...
    """
    try:
        id = re.search(tag + ' "(.+?)"', others).group(1)
    except AttributeError:
        id = ""
    return id
        
def get_novelnum(line):
    """ get novel number such as novel exon, intron, loci from gffcmp.stats """
    return line.rstrip().split("\t")[0].split(" ")[-1].split("/")[0]

def summarize_annotation_update():
    """ summarize annotation update before and after """
    items = ['oldlocinum', 'newlocinum', 'oldlocinum_multitranscript', 'newlocinum_multitranscript', 'oldtranscriptnum', 'newtranscriptnum', 'oldtranscriptnum_multiexon', 'newtranscriptnum_multiexon', 'oldtranscriptnum_perlocus', 'newtranscriptnum_perlocus', 'novellocinum', 'novelexonnum', 'novelintronnum']
    with open("../data/output/annotation.update.summary.tab", 'w') as f:
        f.write("\t" + "\t".join(items) + "\n")
        for dxxx in sharedinfo.ordered_species:
            f.write(dxxx)
            ann = FocalAnnotation(dxxx)
            for item in items:
                f.write("\t" + ann.stats[item])
            f.write("\n")

def summarize_gene_merge():
    """ summarize gene merge number """
    with open("../data/output/gene.merge.summary.tab", 'w') as f:
        for species in sharedinfo.ordered_species:
            ann = FocalAnnotation(species)
            dct_merge = dict()
            dct_merge[1] = 0
            dct_merge[2] = 0
            for newgeneid in ann.new2old_geneid.keys():
                oldgeneid = ann.new2old_geneid[newgeneid]
                print(species,newgeneid,len(oldgeneid),oldgeneid)
                if len(oldgeneid) == 1:
                    dct_merge[1] += 1
                else:
                    dct_merge[2] += 1
            f.write(species + "\t" + str(dct_merge[1]) + "\t" + str(dct_merge[2]) + "\n")

def remove_redundancy_in_dict(dct):
    """ 
       input dct: "a" = {"a1"}, "b" = {"b1", "b2"}
       only one2one key2value were kept in the dict
       output dct: "a" = {"a1"}
    """
    dct_u = dict()
    for i in dct.keys():
        ii = dct[i]
        if len(ii) == 1:
            dct_u[i] = list(ii)[0]
    return(dct_u)

def switch_key_value_in_dict(dct):
    return({y:x for x,y in dct.items()})

class FocalAnnotation:
    """FocalAnnotation object"""
    def __init__(self, species):
        self.species = species
        # self.filename = species + ".SVGpredAdded.v2.gtf"
        self.filename = species + ".v3.gtf"
        self.lines = get_lines("../data/annotation", self.filename)
        self.GdxxxID = self.get_geneid_from_annotation(ver="new")        

        # old annotation (to compare)
        self.oldfilename = species + ".gtf"
        self.oldlines = get_lines("../data/annotation", self.oldfilename)
        self.dxxxFBgn = self.get_geneid_from_annotation(ver="old")

        # compare updated annotation with old (postive v3 vs old-gtf; negative old-gtf vs v3)
        self.statspositivelines = get_lines("../data/annotation", species + ".statspositive")
        self.statsnegativelines = get_lines("../data/annotation", species + ".statsnegative")
        self.stats = self.get_stats()
        self.tmaplines = get_lines("../data/annotation", species + ".v3.gtf.tmap")
        self.dxxxFBgn2GdxxxID, self.GdxxxID2dxxxFBgn = self.get_tmap()

        # there are ortholog one2one (GdxxxID 2 dxxxFBgn 2 dmelFBgn)
        self.filenameolo121 = species + ".concise.one2one"        
        self.linesolo121 = get_lines("../data/ortholog", self.filenameolo121)
        (self.GdxxxID2dmelFBgn, self.dxxxFBgn2dmelFBgn) = self.get_olo121_info()


        # these dmelFBgn lost one2one due to annotation update, and I need to use dxxxFBgn 
        self.lost1_dmelFBgn = set(self.dxxxFBgn2dmelFBgn.values()) - set(self.GdxxxID2dmelFBgn.values())
        self.lost1_dxxxFBgn = set([switch_key_value_in_dict(self.dxxxFBgn2dmelFBgn)[i] for i in self.lost1_dmelFBgn])

        ## there are some newly identified orthologs 121
        self.filenameneo121 = species + ".neo121"
        self.linesneo121 = get_lines("../data/ortholog", self.filenameneo121)
        self.add_neo_to_GdxxxID2dmelFBgn()


        
    def get_geneid_from_annotation(self, ver="new"):
        """ get all geneid as set() """
        lines = self.lines
        if ver == "old":
            lines = self.oldlines

        st = set()
        for line in lines:
            (scaffold, tag, feature, start, end, scoare, strand, dot, others) = line.rstrip().split("\t")
            this_geneid = get_id(others, 'gene_id')
            st.add(this_geneid)
        return(st) 

    def get_new2ref_old_geneid(self):
        """ new2old_geneid
            replace old_geneid (dxxx FBgn) to ref_old_geneid (dmel FBgn)
            by olo121 table
        """
        dct = dict()
        for i in self.new2old_geneid:
            if i in self.olo121.keys():
                ii = self.olo121[i]
                dct[i] = ii
        return(dct)

    def add_neo_to_GdxxxID2dmelFBgn(self):
        """
           add new information (my own orhtolog one2one) to GdxxxID2dmelFBgn
        """
        count = 0
        for line in self.linesneo121:
            (dmelFBgn, GdxxxID) = line.rstrip().split("\t")
            self.GdxxxID2dmelFBgn[GdxxxID] = dmelFBgn
            count += 1

    def get_geneidtable(self):
        """ get geneid, refgeneid, and relationship """
        st_geneid = set()
        st_refgeneid = set()
        dct_geneid2refgeneid = dict()
        dct_refgeneid2geneid = dict()
        for line in self.lines:
            (scaffold, tag, feature, start, end, scoare, strand, dot, others) = line.rstrip().split("\t")
            this_geneid = get_id(others, 'gene_id')
            this_refgeneid = get_id(others, 'ref_gene_id')

            count = 0
            if not this_geneid == "":
                st_geneid.add(this_geneid)
                count += 1
            if not this_refgeneid == "":
                st_refgeneid.add(this_refgeneid)
                count += 1
            if count == 2:
                if not this_geneid in dct_geneid2refgeneid.keys():
                    dct_geneid2refgeneid[this_geneid] = set()
                if not this_refgeneid in dct_refgeneid2geneid.keys():
                    dct_refgeneid2geneid[this_refgeneid] = set()
                
                dct_geneid2refgeneid[this_geneid].add(this_refgeneid)
                dct_refgeneid2geneid[this_refgeneid].add(this_geneid)
        
        return st_geneid, st_refgeneid, dct_geneid2refgeneid, dct_refgeneid2geneid
    
    def get_olo121_info(self):
        """
           get old ortholog from FlyBase/OrthDB
           return many dct
           GdxxxID to dmelFBgn
           dxxxFBgn to dmelFBgn
        """
        dct_GdxxxID2dmelFBgn = dict()
        dct_dxxxFBgn2dmelFBgn = dict()
        for line in self.linesolo121:
            dmelFBgn, dmelSymbol, dxxxFBgn = line.rstrip().split("\t")
            if dxxxFBgn in self.dxxxFBgn2GdxxxID.keys():
                GdxxxID = self.dxxxFBgn2GdxxxID[dxxxFBgn]
                dct_GdxxxID2dmelFBgn[GdxxxID] = dmelFBgn
            dct_dxxxFBgn2dmelFBgn[dxxxFBgn] = dmelFBgn
        return(dct_GdxxxID2dmelFBgn, dct_dxxxFBgn2dmelFBgn)

    def get_121_id(self):
        """ get one to one (geneid to refgeneid)"""

        """ make sure geneid and refgeneid have one to one relationship """
        geneid_one = set()
        for geneid in self.geneid2refgeneid.keys():
            if len(self.geneid2refgeneid[geneid]) == 1:
                geneid_one.add(geneid)

        refgeneid_one = set()
        for refgeneid in self.refgeneid2geneid.keys():
            if len(self.refgeneid2geneid[refgeneid]) == 1:
                refgeneid_one.add(refgeneid)

        dct_geneid2refgeneid121 = dict()
        dct_refgeneid2geneid121 = dict()
        for geneid in geneid_one:
            for refgeneid in self.geneid2refgeneid[geneid]:
                if refgeneid in refgeneid_one:
                    dct_geneid2refgeneid121[geneid] = refgeneid
                    dct_refgeneid2geneid121[refgeneid] = geneid

        """ 
            two manually corrected ortholog
            refgeneid add species to distinguish
        """
        dsx_geneid = sharedinfo.dsx_species2geneid[self.species]
        dsx_refgeneid = sharedinfo.dsx_refgeneid + "_" + self.species
        fru_geneid = sharedinfo.fru_species2geneid[self.species]
        fru_refgeneid = sharedinfo.fru_refgeneid + "_" + self.species

        dct_geneid2refgeneid121[dsx_geneid] = dsx_refgeneid
        dct_refgeneid2geneid121[dsx_refgeneid] = dsx_geneid
        dct_geneid2refgeneid121[fru_geneid] = fru_refgeneid
        dct_refgeneid2geneid121[fru_refgeneid] = fru_geneid

        """
           thousands of newly annotated one to one genes
        """
        for line in self.linesneo121:
            (dmel_refgeneid, dxxx_geneid) = line.rstrip().split("\t")
            dxxx_refgeneid = dmel_refgeneid + "_" + self.species            
            dct_geneid2refgeneid121[dxxx_geneid] = dxxx_refgeneid
            dct_refgeneid2geneid121[dxxx_refgeneid] = dxxx_geneid

        return dct_geneid2refgeneid121, dct_refgeneid2geneid121

    def get_stats(self):
        """ get stats of annotation comparison, such as novel exon, intron, loci """
        dct = dict()
        for line in self.statsnegativelines:
            if "#     Query mRNAs" in line:
                m = re.search('Query mRNAs :   (\d+) in   (\d+) loci  \((\d+) multi-exon', line) 
                dct['oldtranscriptnum'] = m.group(1)
                dct['oldlocinum'] = m.group(2)
                dct['oldtranscriptnum_multiexon'] = m.group(3)
            elif "multi-transcript loci" in line:
                m = re.search('(\d+) multi\-transcript loci, ~(.+?) transcripts', line)
                dct['oldlocinum_multitranscript'] = m.group(1)
                dct['oldtranscriptnum_perlocus'] = m.group(2)

        for line in self.statspositivelines:
            if "#     Query mRNAs" in line:
                m = re.search('Query mRNAs :   (\d+) in   (\d+) loci  \((\d+) multi-exon', line)
                dct['newtranscriptnum'] = m.group(1)
                dct['newlocinum'] = m.group(2)
                dct['newtranscriptnum_multiexon'] = m.group(3)
            elif "multi-transcript loci" in line:
                m = re.search('(\d+) multi\-transcript loci, ~(.+?) transcripts', line)
                dct['newlocinum_multitranscript'] = m.group(1)
                dct['newtranscriptnum_perlocus'] = m.group(2)
            elif "Novel exons" in line:
                dct['novelexonnum'] = get_novelnum(line)
            elif "Novel introns" in line:
                dct['novelintronnum'] = get_novelnum(line)
            elif "Novel loci" in line:
                dct['novellocinum'] = get_novelnum(line)
        return dct

    def get_tmap(self):
        """ 
           get connection of gene ids between gtf and v3.gtf
           only complete one2one relationship will be returned
        """
        old2new_geneid = dict()
        new2old_geneid = dict()
        for line in self.tmaplines:
            (ref_gene_id, ref_id, class_code, qry_gene_id, qry_id, FMI, FPKM, FPKM_conf_lo, FPKM_conf_hi, cov, len, major_iso_id, ref_match_len) = line.rstrip().split("\t")
            if class_code == "=":
                if not ref_gene_id in old2new_geneid.keys():
                    old2new_geneid[ref_gene_id] = set()
                old2new_geneid[ref_gene_id].add(qry_gene_id)
                
                if not qry_gene_id in new2old_geneid.keys():
                    new2old_geneid[qry_gene_id] = set()
                new2old_geneid[qry_gene_id].add(ref_gene_id)

        # make non-redundant version (only one2one were kept in the result)
        old2new_geneid_u = remove_redundancy_in_dict(old2new_geneid)
        new2old_geneid_u = remove_redundancy_in_dict(new2old_geneid)
        return (old2new_geneid_u, new2old_geneid_u)
        
        
if __name__ == '__main__':
    # summarize_annotation_update()
    # summarize_gene_merge()
    ann = FocalAnnotation("dyak")


