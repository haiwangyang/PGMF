#!/usr/bin/env python

"""
Purpose:
    collect coverage of junctional reads in sexed-tissues for certain species
"""
import sharedinfo
import focaljunction
import pandas as pd

def write_junction_summary_to_table(species, locations, genesymbol):
    dcts = list()
    for location in locations:
        dct = focaljunction.get_junction_of_species_by_location(species, location)
        dcts.append(dct)
    mdct = focaljunction.merge_dcts(dcts)
    mpd = pd.DataFrame.from_dict(mdct)

    # reorder columns by fixed order
    mpd = mpd.ix[:, [species + "_" + i for i in sharedinfo.ordered_sexedtissue7]]
    mpd = mpd.reindex(locations)
    mpd.to_csv("../data/output/" + species + "." + genesymbol + ".junc.summary.txt", sep="\t")


if __name__ == '__main__':
    """ check each tra ortholog mannually """
    # focaljunction.get_junction_of_species_by_partiallocation("dmel", "3L:16590")
    # consensus, female-funtional, unisex-dysfunctional1
    # write_junction_summary_to_table("dmel", ["3L:16590295-16590351_-", "3L:16590756-16591003_-", "3L:16590931-16591003_-"], "tra")

    # focaljunction.get_junction_of_species_by_partiallocation("dyak", "3L:17475")
    # consensus, female-funtional, unisex-dysfunctional1
    # write_junction_summary_to_table("dyak", ["3L:17475357-17475426_+", "3L:17474772-17475015_+", "3L:17474772-17474844_+"], "tra")

    # focaljunction.get_junction_of_species_by_partiallocation("dana", "scaffold_13337:13997")
    # consensus, female-funtional, unisex-dysfunctional1, unisex-dysfunctional2
    # write_junction_summary_to_table("dana", ["scaffold_13337:13997573-13997687_+", "scaffold_13337:13996943-13997204_+", "scaffold_13337:13996943-13997054_+", "scaffold_13337:13997139-13997204_+"], "tra")

    # focaljunction.get_junction_of_species_by_partiallocation("dpse", "XR_group8:2283")
    # consensus, female-funtional, unisex-dysfunctional1
    write_junction_summary_to_table("dpse", ["XR_group8:2282529-2282592_-", "XR_group8:2282937-2283184_-", "XR_group8:2283110-2283184_-"], "tra")

    # focaljunction.get_junction_of_species_by_partiallocation("dper", "scaffold_24:772")
    # consensus, female-funtional, unisex-dysfunctional1
    # write_junction_summary_to_table("dper", ["scaffold_24:773156-773218_+", "scaffold_24:772564-772811_+", "scaffold_24:772564-772638_+"], "tra")

    # focaljunction.get_junction_of_species_by_partiallocation("dwil", "scf2_1100000004768:2692")   
    # consensus, female-funtional, unisex-dysfunctional1, unisex-dysfunctional2
    # write_junction_summary_to_table("dwil", ["scf2_1100000004768:2692832-2692893_+", "scf2_1100000004768:2692173-2692397_+", "scf2_1100000004768:2692173-2692248_+", "scf2_1100000004768:2692308-2692397_+"], "tra")

    # focaljunction.get_junction_of_species_by_partiallocation("dmoj", "scaffold_6654:1870")
    # consensus, female-funtional, unisex-dysfunctional1, unisex-dysfunctional2
    # write_junction_summary_to_table("dmoj", ["scaffold_6654:1870972-1871031_+", "scaffold_6654:1870338-1870585_+", "scaffold_6654:1870338-1870407_+", "scaffold_6654:1870518-1870585_+"], "tra")

    # focaljunction.get_junction_of_species_by_partiallocation("dvir", "scaffold_13049:1963")
    # consensus, female-funtional, unisex-dysfunctional1
    write_junction_summary_to_table("dvir", ["scaffold_13049:1962579-1962641_-", "scaffold_13049:1963022-1963278_-", "scaffold_13049:1963210-1963278_-"], "tra")
 
    # focaljunction.get_junction_of_species_by_partiallocation("dgri", "scaffold_15110:3027")
    # consensus, female-funtional, unisex-dysfunctional1
    write_junction_summary_to_table("dgri", ["scaffold_15110:3026921-3027142_-","scaffold_15110:3027520-3027803_-", "scaffold_15110:3027724-3027803_-"], "tra")
    # p3 funtionaltra.py | g 3027520
