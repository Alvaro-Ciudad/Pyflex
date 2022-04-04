#!/usr/bin/env python

# IMPORTS
import argparse
import subprocess
from getopt import getopt
from Bio import PDB
from sympy import sequence
import os
import sys

sequencee = "MALLRDVSLQDPRDRFELLQRVGAGTYGDVYKARDTVTSELAAVKIVKLDPGDDISSLQQEITILRECRHPNVVAYIGSYLRNDRLWICMEFCGGGSLQEIYHATGPLEERQIAYVCREALKGLHHLHSQGKIHRDIKGANLLLTLQGDVKLADFGVSGELTASVAKRRSFIGTPYWALMLMSKSSFQPPKLRDKTRWTQNFHHFLKLALTKNPKKRPTAEKLLQHPFTTQQLPRALLTQLLDKASDPHLGTPSPEDCELETYDMFPDTIHSRGQHGPAERTPSEIQFHQVKFGAPRRKETDPLNEP"


# PSI-BLAST
def get_psiblast_prot(sequence, database):
    """This function returns the result of the psi-blast REST API"""
    os.chdir('./blastscript/')
    p = subprocess.Popen(['python', 'psiblast.py', '--email', 'pyflex@protonmail.com',
                          '--sequence', sequence, '--database', database], stdout=subprocess.PIPE)
    print("Psi-blast is running...")
    i = ""
    for line in p.stdout:
        i = line
    p.wait()
    print(p.returncode)
    os.chdir('../')  # IMPORTANT

    return i


def decode_job_id():
    """This function returns the job ID as a string, because the result of the previously subprocess is always a bytes object"""
    bytedata = get_psiblast_prot()
    job_id = bytedata.decode('UTF-8')
    job_id = job_id.split(" ")
    job_id = job_id[3].split(".")
    return job_id[0]


p = get_psiblast_prot(sequencee, "uniprotkb")
d = decode_job_id()


print(p)
print(d)
