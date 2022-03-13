#!/usr/bin/env python


# IMPORTS
import argparse
import subprocess
from getopt import getopt
from Bio import PDB
from prody import *
from matplotlib.pylab import *

# imports from psiblast.py REST API
import psiblast as pb






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

options = parser.parse_args()


print(options.infile)


# INPUT
def get_input_prot():
    """This function returns the input protein depending on the introduced way to the program."""


    protein = "VAYIGSYLRNDRLWICMEFCGGGSLQEIYHATGPLEERQ"
    return protein

# PSI-BLAST
def get_psiblast_prot():
    """This function returns the result of the psi-blast REST API"""
    # command line
    cmd_l = "python psiblast.py --email pyflex@protonmail.com --sequence MALLRDVSLQDPRDRFELLQRVGAGTYGDVYKARDTVTSELAAVKIVKLDPGDDISSLQQEITILRECRHPNVVAYIGSYLRNDRLWICMEFCGGGSLQEIYHATGPLEERQIAYVCREALKGLHHLHSQGKIHRDIKGANLLLTLQGDVKLADFGVSGELTASVAKRRSFIGTPYWALMLMSKSSFQPPKLRDKTRWTQNFHHFLKLALTKNPKKRPTAEKLLQHPFTTQQLPRALLTQLLDKASDPHLGTPSPEDCELETYDMFPDTIHSRGQHGPAERTPSEIQFHQVKFGAPRRKETDPLNEP --database uniprotkb"
    p = subprocess.Popen(cmd_l, stdout=subprocess.PIPE)
    for line in p.e:
        print(line)
    p.wait()
    print(p.returncode)


    

#seq = get_input_prot()
psi_blast_result = get_psiblast_prot()



# DB ALPHAFOLD

# HOMOLOGY

# DISTANCES



