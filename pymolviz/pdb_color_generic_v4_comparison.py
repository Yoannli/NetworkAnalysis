#!/usr/bin/env python
#
"""
diverging palette explained:
If the value 0 is meaningful, then 
"""


'''pdb_color_generic_v2.py adapted from Warren Francis' https://github.com/wrf/pdbcolor

    options are:
    -s sequence names to get the chain from DBREF records, otherwise use default chain
    -c column, column index of desired data, starting from 0
    -d delimiter, default is tab, e.g. change to "," for csv
    -g group name, changes the group name in the PyMOL selection
    -l color schemes, schemes include:
      sequential - from gray to intense color (e.g. -l green )
        red, yellow, blue, green
      diverging - with either light or dark in the middle ( -l div2b )
      so values close to 0 will blend in with the background color
        div1w - dark brown to white to green, for white bg display
        div1b - light brown to black to light green, for black bg display
        div2w - brick red to white to dark blue, for white bg
        div2b - pink to black to sky blue, for black bg
    --exclude-common do not print the bin for the most common group
    -O specify which group is common, otherwise assumed to be either 0, 4, or 8
    -f, residue number(s) to highlight with a different color

    data -i can be any text file, so long as columns can be split with -d
    generally, it is best to format in something like this:
residue score   chain
1       88.5    A

    chain information is optional, and the residue column can be specified
    with the option --site-column, starting from 0

    if multiple proteins are in the PDB file, then it is better to have
    the chain specified in the data file. Otherwise the script can be run
    multiple times, each time specifying a different chain by --default-chain
'''

import sys
import argparse
import itertools
import numpy # for arange,around,floor,log10,median
import re
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import pickle # to save the palette only


def read_generic_data(datafile, delimiter, scorecolumn, pdbcolumn=0, sitecolumn=1, chaincolumn=2, defaultchain="A"):
    '''read generic data file, and return a dict key is site number and value is score'''
    rawscore_dict = dict() 
    linecounter = 0
    foundscores = 0
    columnmax = 0

    sys.stderr.write("\n\n=========================================================\n")
    sys.stderr.write("# Reading scores from column {} in {}, separated by {}\n".format( scorecolumn, datafile, delimiter ) )

    with open(datafile, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line and line[0]!="#":
                linecounter += 1
                lsplits = line.split(delimiter)
                if len(lsplits) > columnmax:
                    columnmax = len(lsplits)
                sitecol = lsplits[sitecolumn] # raw, might be string
                try: # check if residue number can be turned into integer
                    sitecol = re.sub("[^0-9]", "", sitecol) # remove non-numeric characters
                    site = int(float(sitecol))
                except ValueError: # no value, row is empty in this column or is header
                    continue
                try: # check if score can be turned into integer
                    # score = re.sub("[^0-9.]", "", lsplits[scorecolumn]) # remove non-numeric characters
                    score = lsplits[scorecolumn]
                    score = float(score)
                except ValueError: # no value, row is empty in this column or is header
                    continue
                if type(chaincolumn) is int:
                    chain = lsplits[int(chaincolumn)]
                else: # use default chain
                    chain = defaultchain
                pdb = lsplits[pdbcolumn]
                if pdb not in rawscore_dict:
                    rawscore_dict[pdb] = defaultdict(dict) # key is chain, value is dict of site and raw score
                rawscore_dict[pdb][chain][site] = float(score)
                foundscores += 1		
    sys.stderr.write("# Counted {} lines with {} scores\n".format( linecounter, foundscores ) )
    if foundscores==0:
        if columnmax < 2:
            sys.stderr.write("# WARNING: No scores found, lines contained only 1 column, check -d\n")
    return rawscore_dict

def get_chains_only(defaultchain, seqidlist, pdbfile):
    '''read PDB file and return two dicts, one where key is chain and value is sequence ID, other where key is the chain and value is integer of the DBREF offset'''
    keepchains = {} # dict where key is chain and value is seqid, though value is not used
    refoffsets = {} # key is chain, value is integer offset from DB seq
    sys.stderr.write("# Reading chain from PDB {}\n".format(pdbfile) )
    for line in open(pdbfile,'r'):
        record = line[0:6].strip()
        # get relevant chains that match the sequence, in case of hetero multimers
        if record=="DBREF":
            sys.stderr.write(f"\n### found DBREF record {seqidlist}\n\n")
            defaultchain = False
            proteinid = line[42:56].strip()
            for seqid in seqidlist:
                if seqid.find(proteinid)>-1:
                    chaintarget = line[12]
                    chainstart = int(line[14:18].strip())
                    dbstart = int(line[55:60].strip())
                    chainoffset = dbstart - chainstart
                    sys.stderr.write("### keeping chain {} for sequence {} with offset {}\n".format( chaintarget, proteinid, chainoffset ) )
                    keepchains[chaintarget] = proteinid
                    refoffsets[chaintarget] = chainoffset
    if defaultchain: # meaning nothing was found, use default and single sequence
        if seqidlist: # all default chains are assumed to use only the first sequence
            keepchains[defaultchain] = seqidlist[0]
        else: # value is not called, but just to indicate that -s is or used or not
            keepchains[defaultchain] = "UNKNOWN"
        refoffsets[defaultchain] = 0
        sys.stderr.write("### using default chain {}\n".format( defaultchain ) )
    return keepchains, refoffsets


def create_diverging_palette(data, meaninful_zero=True, num_bins=11, palette="BrBG"):
    # create a diverging palette based on the data. Assumes data spans 0. 
    # if meaninful_zero is True, then 0 is the centre of the data/palette. If False the centre is the mean of the data

    tmp_num_bins = num_bins - 1  # for calculating low_set and high_set (actual number of bins will be original num_bins)
    min_val, max_val = min(data), max(data)
    val_range = max_val - min_val

    # create bins with equal range
    bin_range = val_range / tmp_num_bins

    if meaninful_zero: # for example functional scores
        low_mid_bin = 0 - bin_range
        high_mid_bin = 0 + bin_range

        binvalues = [low_mid_bin, 0, high_mid_bin]

        lowest = low_mid_bin
        while lowest > min_val:
            lowest -= bin_range
            binvalues.insert(0, lowest)
        highest = high_mid_bin
        while highest < max_val:
            highest += bin_range
            binvalues.append(highest)

        n_less_than_zero = len([x for x in binvalues if x < 0])  # number of bins in negative range
        n_more_than_zero = len([x for x in binvalues if x > 0])  # number of bins in positive range
        if abs(min_val) > abs(max_val):
            colour_maps = np.vstack((
                plt.get_cmap(palette)(np.linspace(0, 0.5, n_less_than_zero + 1))[:-1],
                plt.get_cmap(palette)(np.linspace(0.5, 1, n_less_than_zero + 1))[1:1 + n_more_than_zero]
            ))
        else:
            colour_maps = np.vstack((
                plt.get_cmap(palette)(np.linspace(0, 0.5, n_more_than_zero + 1))[-1 - n_less_than_zero:-1],
                plt.get_cmap(palette)(np.linspace(0.5, 1, n_more_than_zero + 1))[1:]
            ))

        colour_maps = [x[:-1] for x in colour_maps] # remove alpha channel
        colour_maps = [[round(x, 3) for x in y] for y in colour_maps]
        binvalues = [round(x, 2) for x in binvalues]
        binvalues[0] -= 0.01 # to make sure it includes the lowest value
        binvalues[-1] += 0.01 # to make sure it includes the highest value
        return colour_maps, binvalues

    else: # for example network scores, where 0 does not mean anything
        binvalues = np.linspace(min_val, max_val, num_bins+1)
        binvalues = [round(x, 2) for x in binvalues]
        binvalues[0] -= 0.01 # to make sure it includes the lowest value
        binvalues[-1] += 0.01 # to make sure it includes the highest value
        colour_maps = plt.get_cmap(palette)(np.linspace(0, 1, num_bins))
        colour_maps = [x[:-1] for x in colour_maps] # remove alpha channel
        colour_maps = [[round(x, 3) for x in y] for y in colour_maps]
        return colour_maps, binvalues




def make_output_script(wayout, scoredict, groupname, palette="BrBG", reverse_colors=False, highlight_residue=None, meaninful_zero=True):

    try:
        plt.get_cmap(palette)
    except:
        sys.stderr.write(f"Invalid palette name: {palette}. Switching to default palette: BrBG.")
        palette = "BrBG"

    default_color = plt.get_cmap(palette)(0.5)[:-1]  # default color for residues not in any group

    # print(scoredict)
    
    # all_scores = list( itertools.chain( list(rd.values()) for rd in scoredict.values() ) )[0]

    all_scores = [score for pdb in scoredict for chain in scoredict[pdb] for site in scoredict[pdb][chain] for score in [scoredict[pdb][chain][site]]]
    # print(all_scores)
    sys.stderr.write("# Generating list of bins from data\n")
    lowest_val = min(all_scores)
    highest_val = max(all_scores)
    val_range = highest_val - lowest_val
    median_val = numpy.median(all_scores)
    sys.stderr.write("# data range from {:.2f} to {:.2f}, diff of {:.2f}, median of {:.2f}\n".format(lowest_val, highest_val, val_range, median_val) )

    NUM_BINS = 11
    targetcolors, binvalues = create_diverging_palette(all_scores, meaninful_zero=meaninful_zero, num_bins=NUM_BINS, palette=palette)
    if reverse_colors:
        targetcolors.reverse()
    # save targetcolors and binvalues to a file
    with open("palette.pkl", "wb") as f:
        pickle.dump((targetcolors, binvalues), f)


    sys.stderr.write("# Generating PyMOL script for color scheme {} with bins of:\n{}\n".format( palette, binvalues ) )

    # set color for all objects that are not part of the target chains
    wayout.write("set_color colordefault, [{}]\n".format( ",".join(map(str,default_color)) ) )
    wayout.write("color colordefault, all\n")

    alphabets = 'abcdefghijklmnopqrstuvwxyz'
    colour_levels = [x for x in alphabets[:len(binvalues)-1]]

    for i,rgb in enumerate(targetcolors):
        # note: had to correct this since bins for -1.5 and -1.0 were being assigned the same color name since it uses int 
        # colorname = "{}{:02d}".format( basecolor, int(binvalues[i]*binname_correction) ) # previous code
        colorname = "{}_{}".format( palette, colour_levels[i] )
        wayout.write("set_color {}, [{}]\n".format( colorname, ",".join(map(str,rgb)) ) )
    if highlight_residue is not None:
        highlight_color = [1.0, 0.0, 0.0]  # You can change this to the desired RGB highlight color
        wayout.write(f"set_color highlight, [{','.join(map(str,highlight_color))}]\n")
    
    # print("\n===================\n")
    # all_scoregroups = {} # for each pdb_chain (key), store the residues in each group (a )
    # make commands for each pdb and chain
    colourgroups = defaultdict(list)
    for pdb in scoredict.keys():
        for chain in scoredict[pdb].keys(): # keys are chain letters, values are seq IDs
            # print(chain)
            scoregroups = defaultdict(list) # key is percent group, value is list of residues
            # for each residue, assign to a bin
            # print(pdb, chain)
            for residue in scoredict[pdb][chain].keys():
                residuescore = scoredict[pdb][chain].get(residue,0.00)
                for i,value in enumerate(binvalues[:-1]):
                    upper = binvalues[i+1]
                    if residuescore < upper:
                        # scoregroups[value].append(residue - chainoffset)
                        scoregroups[value].append(residue)  # chainoffset removed
                        break
            # all_scoregroups[f"{pdb}_{chain}"] = scoregroups
        

            # assign whole chain to lowest color, then build up
            wayout.write("color colordefault, {} & chain {}\n".format( pdb, chain ) )
            for i,value in enumerate(binvalues[:-1]):
                binname = f"{colour_levels[i]}_{groupname}_{i+1}_{chain}"
                resilist = list(map(str,scoregroups[value]))
                if resilist: # do not print empty groups
                    colourgroups[binname].append( (pdb, chain, resilist) )

                    # binresidues = ",".join(resilist)
                    # wayout.write("select {}, ({} & chain {} & resi {})\n".format( binname, pdb, chain, binresidues ) )
                    # wayout.write("color {}_{}, {}\n".format( palette, colour_levels[i], binname ) )
        
            # # Highlight specific residue with a different (currently red) color
            # if highlight_residue is not None:
            #     res_to_highlight = ",".join(map(str, highlight_residue))
            #     wayout.write("select highl, (chain {} & resi {})\n".format( chain, res_to_highlight ) )
            #     wayout.write("color highlight, highl\n")
    # print(colourgroups)
    for binname, groupresidues in colourgroups.items():
        pdbchains = []
        for pdb, chain, reslist in groupresidues:
            pdbchains.append(f"({pdb} & chain {chain} & resi {','.join(reslist)})")
        wayout.write("select {}, ({})\n".format( binname, " or ".join(pdbchains) ) )
        wayout.write("color {}_{}, {}\n".format( palette, binname.split("_")[0], binname ) )

    wayout.write("deselect\n") 
    # no return

def main(argv, wayout):
    if not len(argv):
        argv.append('-h')
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument("-c","--data-column", default=1, type=int, help="index of data column, starting from 0 [1]")
    parser.add_argument("-d","--delimiter", default="\t", help="delimiter for file, default is tab")
    parser.add_argument("-g","--group-name", default="grp", help="name for groups, default is grp, appears as 100_grp_9_A")
    parser.add_argument("-i1","--input-file", help="tabular or csv file of first site data, with sites in the first column", required=True)
    parser.add_argument("-l","--base-color", default="red", help="color gradient, default is red, options are: [red,yellow,blue,green,div1w,div1b,div2w,div2b]")
    # parser.add_argument("-s","--sequence", nargs="*", help="sequence ID for PDB, give multiple names if data is available in the input file")
    parser.add_argument("--default-chain", default="A", help="default letter of chain [A], if DBREF for the sequence cannot be found in PDB")
    parser.add_argument("-x","--exclude-common-group", action="store_true", help="exclude common group, for cases where there are a large number of score-0 residues")
    parser.add_argument("-O","--default-chain-override", type=int, help="index to color all residues by default (from 0 to 8) for lowest score group, otherwise determined automatically")
    parser.add_argument("-r","--reverse-colors", action="store_true", help="if used, reverse colors for negative-value datasets")
    parser.add_argument("-f", "--highlight-residue", nargs="*", type=int, help="residue number to highlight with a different color")
    parser.add_argument("-z","--meaningful-zero", action="store_false", help="if used, middle color is the mean of the data, otherwise it is 0 (default)")
    parser.add_argument("--chain-column", type=int, help="column containing chain ID, default is None, will use chain A")
    parser.add_argument("--site-column", default=1, type=int, help="index of site column, starting from 0 [0]")
    parser.add_argument("--zero-override", default=0.0, type=float, help="middle index if data spans negative to positive, default is 0 as the middle color")
    parser.add_argument("--pdb-column", default=0, type=int, help="column index of PDB ID, default is 0")
    
    args = parser.parse_args(argv)
    # read generic format data
    datadict = read_generic_data(args.input_file, args.delimiter, args.data_column, args.pdb_column, args.site_column, args.chain_column, args.default_chain)

    # make PyMOL script with color commands-
    # refchains, refoffsets = get_chains_only(args.default_chain, args.sequence, args.pdb1) # NOTE WE DONT NEED THIS ANYMORE because we already assumed we filtered for relevant chain
    make_output_script(wayout, datadict, args.group_name, args.base_color, args.reverse_colors, args.highlight_residue, args.meaningful_zero)
    

if __name__ == "__main__":
    main(sys.argv[1:], sys.stdout)
