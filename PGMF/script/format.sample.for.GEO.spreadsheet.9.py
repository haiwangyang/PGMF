from sharedinfo import get_M2N


ts2tissue = get_M2N("../data/tissue/ts2tissue2FBID.tab", 0, 1)
ts2FBID = get_M2N("../data/tissue/ts2tissue2FBID.tab", 0, 2)
sx2sex = get_M2N("../data/sex/sx2sex.tab", 0, 1)
sx2FBID = get_M2N("../data/sex/sx2sex.tab", 0, 2)
sp2species = get_M2N("../data/species/species2dxxx.txt", 0, 2)
sp2ID = get_M2N("../data/species/species2dxxx.txt", 0, 3)
sp2FBsp = get_M2N("../data/species/species2dxxx.txt", 0, 4)

GSM2sample = get_M2N("/Users/yangh13/PGMF/PGMF/data/GSM/GSE99574/GSM2sample", 0, 1)
GSM2SRX = get_M2N("/Users/yangh13/PGMF/PGMF/data/GSM/GSE99574/GSM2SRX.txt", 0, 1)
SRX2SRR = get_M2N("/Users/yangh13/PGMF/PGMF/data/GSM/GSE99574/SRX2SRR.txt", 0, 1)


for GSM in sorted(GSM2sample.keys()):
    sample = GSM2sample[GSM]
    SRX = GSM2SRX[GSM]
    SRR = SRX2SRR[SRX]
    if "ERCC" in sample:
        print("\t".join([GSM, sample]))
    elif "leftover" in sample:
        print("\t".join([GSM, sample]))
    else:
        sp, ts, sx, rep = sample.split("_")
        ts = ts.lower()
        tissue = ts2tissue[ts]
        tsID = ts2FBID[ts]
        sex = sx2sex[sx]
        sxID = sx2FBID[sx]
        species = sp2species[sp]
        spID = sp2ID[sp]
        spFBsp = sp2FBsp[sp]
        print("\t".join([GSM, sample]))



