from sharedinfo import get_M2N

ts2tissue = get_M2N("../data/tissue/ts2tissue2FBID.tab", 0, 1)
ts2FBID = get_M2N("../data/tissue/ts2tissue2FBID.tab", 0, 2)
sx2sex = get_M2N("../data/sex/sx2sex.tab", 0, 1)
sx2FBID = get_M2N("../data/sex/sx2sex.tab", 0, 2)
sp2species = get_M2N("../data/species/species2dxxx.txt", 0, 2)
sp2ID = get_M2N("../data/species/species2dxxx.txt", 0, 3)
sp2FBsp = get_M2N("../data/species/species2dxxx.txt", 0, 4)

GSM2SRX = get_M2N("../data/GSM/GSE.plus.txt", 0, 3)
GSM2SRR = get_M2N("../data/GSM/GSE.plus.txt", 0, 4)
GSM2htseq = get_M2N("../data/GSM/GSE.plus.txt", 0, 5)

GSM2SAMN = get_M2N("../data/GSM/GSE.plus.txt", 0, 6)

with open("../data/sample/sample2GSM.tab", "r") as f:
    for line in f.readlines():
        sample, GSM = line.rstrip().split("\t")
        SRX = GSM2SRX[GSM]
        SRR = GSM2SRR[GSM]
        SAMN = GSM2SAMN[GSM]
        htseq = "-"
        if GSM in GSM2htseq.keys():
            htseq = GSM2htseq[GSM]
        sp, sx, ts, rep = sample.split("_")
        species = sp2species[sp]
        speciesID = sp2ID[sp]
        Dxxx = sp.replace("w1118", "dmel").replace("oreR", "dmel").replace("dgriG1","dgri").title()
        speciesFBsp = sp2FBsp[sp]
        sex = sx2sex[sx]
        sexFBID = sx2FBID[sx]
        tissue = ts2tissue[ts]
        tissueFBID = ts2FBID[ts]
        replicate = rep.replace("R", "")
        annVer = "2017_03"
        if sp == "dgriG1":
            annVer = "2016_05"
        if ts == "re":
            if sp == "w1118" or sp == "oreR" or sp == "dgriG1":
                tissueFBID = "FBbt:00004857 - FBbt:00004858"
                tissue = "internal " + tissue
        if sp == "w1118":
            Dxxx = sp
        elif sp == "oreR":
            Dxxx = sp

        #print(Dxxx + " " + sex + " " + tissue + " replicate " + replicate )
        #print("NCBI BioSample:" + SAMN)
        #print("FlyBase CV:FBbt:00003004, " + sexFBID + ", and " + tissueFBID)
        #print("biological sample (" + species + " " + sex + ")")
        #print("large-scale dataset (" + species + " " + sex + ")")
        #print(Dxxx + " " + sex + " " + tissue + " replicate " + replicate + " gene-level read counts based on STAR alignments and FlyBase " + annVer + " annotation")
        #print("NCBI GEO " + GSM + ":" + htseq)

        
        #print(Dxxx + " " + sex + " " + tissue + " replicate " + replicate + " gene-level read counts based on HiSAT2 alignments and YO annotation")
        #print("NCBI GEO " + GSM + ":" + htseq.replace("txt", "HiSAT2.YO.txt"))

        #print(Dxxx + " " + sex + " " + tissue + " replicate " + replicate + " gene-level read counts based on HiSAT2 alignments and FlyBase " + annVer + " annotation")         
        #print("NCBI GEO " + GSM + ":" + htseq.replace("txt", "HiSAT2.FB.txt"))

        #print(Dxxx + " " + sex + " " + tissue + " replicate " + replicate + " junction reads")

for sp in ['dana', 'dgri', 'dmel', 'dmoj', 'dper', 'dpse', 'dvir', 'dwil', 'dyak']:
    species = sp2species[sp]
    speciesID = sp2ID[sp]
    speciesID = sp2ID[sp]
    Dxxx = sp.title()
    speciesFBsp = sp2FBsp[sp]
    annVer = "2017_03"
    if sp == "dgri":
        annVer = "2016_05"
    #print(Dxxx + " gene-level normalized read counts based on HiSAT2 alignments and YO annotation")
    #print(Dxxx + " gene-level normalized read counts based on HiSAT2 alignments and FlyBase " + annVer + " annotation")
    #print(Dxxx + " transcript-level TPM based on HiSAT2 alignments and YO annotation")
    #print(Dxxx + " transcript-level TPM based on HiSAT2 alignments and FlyBase " + annVer + " annotation")
    #print(Dxxx + " 1:1 Dmel ortholog table based on YO annotation")
    #print(Dxxx + " YO annotation in gtf format based on HiSAT2 alignments")
    print(Dxxx + " YO annotation in gff3 format based on HiSAT2 alignments")
    for dir in ['forward', 'reverse']:
        #print(Dxxx + " species-level bigWig track in " + dir + " strand")    
        pass

with open("../data/sample/sample_no_rep.list", "r") as f:
    for line in f.readlines():
        sample= line.rstrip().split("\t")[0]
        sp, sx, ts = sample.split("_")
        species = sp2species[sp]
        speciesID = sp2ID[sp]
        sex = sx2sex[sx]
        sexFBID = sx2FBID[sx]
        tissue = ts2tissue[ts]
        tissueFBID = ts2FBID[ts]
        if ts == "re":
            if sp == "w1118" or sp == "oreR" or sp == "dgriG1":
                tissueFBID = "FBbt:00004857 - FBbt:00004858 - FBbt:00004829"
                tissue = "internal " + tissue
        for dir in ['forward', 'reverse']:
            # print("Expression track (bigWig); strain (" + speciesID + "); adult (FBbt:00003004); " + sex + " (" + sexFBID + "); " + tissue + " (" + tissueFBID + "); strand (" + dir + ")")
            # print("https://doi.org/10.6084/m9.figshare.6041888.v1:" + sample + "." + dir + ".bw")
            # print("sequence-based reagent (" + species + " " + sex + ")")
            pass

for sp in ['dana', 'dgri', 'dmel', 'dmoj', 'dper', 'dpse', 'dvir', 'dwil', 'dyak']:
    species = sp2species[sp]
    speciesID = sp2ID[sp]
    # print("sequence-based reagent (" + species + ")")
    # print("Gene-level normalized read counts (txt); strain (" + speciesID + "); adult (FBbt:00003004)")
    # print("Ortholog and ID conversion (txt); strain (" + speciesID + "); adult (FBbt:00003004)")
    # print("Updated annotatiohn (gtf); strain (" + speciesID + "); adult (FBbt:00003004)")
    # print("Updated annotatiohn (gff3); strain (" + speciesID + "); adult (FBbt:00003004)")
    # print("https://doi.org/10.6084/m9.figshare.6042005.v1:" + sp.title() + ".YO.gtf")
    # print("https://doi.org/10.6084/m9.figshare.6042005.v1:" + sp.title() + ".YO.gff3")
    # print("large-scale dataset (" + species + ")")
    # print("GEO Illumina RNA-seq dataset; strain (" + speciesID + "); adult (FBbt:00003004)")
    if sp != "dgri":
        pass
        #print("GSE99574:" + sp + "." + ".gene_level_nrc.YO.tar")
        #print("GSE99574:" + sp + "." + ".gene_level_nrc.FB.tar")
        #print("GSE99574:" + sp + "." + ".transcript_level_TPM.YO.tar")
        #print("GSE99574:" + sp + "." + ".transcript_level_TPM.FB.tar")
        #print("GSE99574:" + sp + ".ortholog_id_conversion.tar")
        #print("GSE99574:" + sp + ".YO.gff3.tar")
    else:
        pass
        #print("GSE80124:" + sp + "." + ".gene_level_nrc.YO.tar")
        #print("GSE80124:" + sp + "." + ".gene_level_nrc.FB.tar")
        #print("GSE80124:" + sp + "." + ".transcript_level_TPM.YO.tar")
        #print("GSE80124:" + sp + "." + ".transcript_level_TPM.FB.tar")
        #print("GSE80124:" + sp + ".ortholog_id_conversion.tar")
        #print("GSE80124:" + sp + ".YO.gff3.tar")
    #print("; ".join(["gene-level normalized read counts based on YO annotation(txt,tar)", species+"(" + speciesFBsp + ")", "strain(" + speciesID + ")", "adult(FBbt:00003004)"]))
    #print("; ".join(["gene-level normalized read counts based on FlyBase annotation(txt,tar)", species+"(" + speciesFBsp + ")", "strain(" + speciesID + ")", "adult(FBbt:00003004)"]))
    #print("; ".join(["transcript-level TPM based on YO annotation(txt,tar)", species+"(" + speciesFBsp + ")", "strain(" + speciesID + ")", "adult(FBbt:00003004)"]))
    #print("; ".join(["transcript-level TPM based on FlyBase annotation(txt,tar)", species+"(" + speciesFBsp + ")", "strain(" + speciesID + ")", "adult(FBbt:00003004)"]))
    #print("; ".join(["1:1 ortholog table(txt,tar)", species+"(" + speciesFBsp + ")", "strain(" + speciesID + ")", "adult(FBbt:00003004)"]))
    #print("; ".join(["YO annotation(gtf,tar)", species+"(" + speciesFBsp + ")", "strain(" + speciesID + ")", "adult(FBbt:00003004)"]))
    #print("; ".join(["YO annotation(gff3,tar)", species+"(" + speciesFBsp + ")", "strain(" + speciesID + ")", "adult(FBbt:00003004)"]))
    #print("; ".join(["GEO Illumina RNA-seq dataset", species+"(" + speciesFBsp + ")", "strain(" + speciesID + ")", "adult(FBbt:00003004)"]))
    for dir in ['forward', 'reverse']:
        pass
        #print("sequence-based reagent (" + species + ")")
        #print("; ".join(["bigWig track(tar)", species+"(" + speciesFBsp + ")", "strain(" + speciesID + ")", "adult(FBbt:00003004)", "strand(" + dir + ")"]))
        if sp != "dgri":
            pass
            #print("GSE99574:" + sp + "." + dir + ".bigWig.tar")
        else:
            pass
            #print("GSE80124:" + sp + "." + dir + ".bigWig.tar")
