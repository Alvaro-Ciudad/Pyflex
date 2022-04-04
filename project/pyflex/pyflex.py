#!/usr/bin/env python


# IMPORTS
import argparse
import subprocess
from getopt import getopt
from Bio import PDB
from prody import *
from matplotlib.pylab import *

from putupsiblast import get_psiblast_prot


# COMMAND LINE ARGUMENTS
parser = argparse.ArgumentParser(description="This program does blablablabla")

requiredNamed = parser.add_argument_group("Required arguments")

requiredNamed.add_argument('-i', '--input',
                    dest = "infile",
                    action = "store",
                    default = None,
                    required = True,
                    help = "Input FASTA file/Uniprot Id")

requiredNamed.add_argument('-f', '--format',
                            required = True,
                            help = "Input format")

parser.add_argument('-o', '--output',
                    dest = "outfile",
                    action = "store",
                    default = "output.pdb",
                    help = "Output file")

parser.add_argument('-g', '--graph',
                    dest = "outgraph",
                    action = "store",
                    default = None,
                    help = "Graphical result")

(options, args) = parser.parse_args()

print(options.infile)







# DB ALPHAFOLD

# HOMOLOGY

# DISTANCES
