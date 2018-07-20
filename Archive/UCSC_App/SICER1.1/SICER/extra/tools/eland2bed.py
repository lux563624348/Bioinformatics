#!/usr/bin/env python
# Copyright (c) 2007 NHLBI, NIH
# Authors: Dustin E Schones and Keji Zhao
#
# This software is distributable under the terms of the GNU
# General Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or otherwise
# using this module constitutes acceptance of the terms of this License.
#
# Disclaimer
# 
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# schonesde@mail.nih.gov).
# 
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser



""" In this version of this program the coords for each match type are
stored in arrays - if there are duplicates, they are inlcluded, but
they are not noted as duplicates in any way -- this was taking too
long

Eventually this should all be updated to use the BED classes.

"""


U0_coords = [];
U1_coords = [];
U2_coords = [];



def getThisCoord(sline, match_id):
    seqlen = len(sline[1]);
    strand = '';
    score = 0
    coord = re.sub(".fa", "", sline[6]);
    start = atoi(sline[7]) - 1;  ## -1 to deal with Eland starting at 1 and UCSC starting at 0
    stop = start + seqlen;

    coord += "\t" + str(start) + "\t" + str(stop) + "\t" + match_id + "\t" + str(score) + "\t";
    if re.match("F", sline[8]):
        coord += "+"
    elif re.match("R", sline[8]):
        coord += "-";
    return coord;


def getCoords(file):
    infile = open(file, 'r');
    for line in infile:
        line = line.strip();
        sline = line.split();
        if re.match("U0", sline[2]):
            coord = getThisCoord(sline, "U0");
            U0_coords.append(coord);
        elif re.match("U1", sline[2]):
            coord = getThisCoord(sline, "U1");
            U1_coords.append(coord);
        elif re.match("U2", sline[2]):
            coord = getThisCoord(sline, "U2");
            U2_coords.append(coord);


def printBED(out):
    outfile = open(out, 'w');
    for c in U0_coords:
        outline = c + "\n";
        outfile.write(outline);
    for c in U1_coords:
        outline = c + "\n";
        outfile.write(outline);
    for c in U2_coords:
        outline = c + "\n";
        outfile.write(outline);


def main(argv):
    parser = OptionParser()
    parser.add_option("-r", "--result_file", action="store", type="string",
                      dest="resultfile", help="eland result file", metavar="<file>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output file name", metavar="<file>")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 4:
        parser.print_help()
        sys.exit(1)

    getCoords(opt.resultfile);
    printBED(opt.outfile);
    

if __name__ == "__main__":
    main(sys.argv)


        
