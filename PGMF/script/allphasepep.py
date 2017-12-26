import focalgene

def get_species2geneid(genesymbol):
    """ get a dict of species to geneid for certain genesymbol """
    dct = {}
    with open("../data/ortholog/" + genesymbol + ".species2geneid", 'r') as f:
        for line in f.readlines():
            species, geneid = line.rstrip().split("\t")
            dct[species] = geneid
    return dct

def print_out_all_possible_pep(species, genesymbol, geneid):
    fg = focalgene.FocalGene(species, geneid, "M")
    for transid in fg.isoform.keys():
        isoformdanseq = fg.isoform[transid]['isoformdnaseq']
        for i in range(0,3):
            print(">" + species + "_" + transid + "_" +  "phase" + str(i) + " " + genesymbol + " " + geneid + "\n" + focalgene.get_pep_and_leftover_from_dna(isoformdanseq, i)[0])

if __name__ == '__main__':
    for genesymbol in ('tra',):
        dsx_dct = get_species2geneid(genesymbol)
        for species, geneid in dsx_dct.items():
            print_out_all_possible_pep(species, genesymbol, geneid)

