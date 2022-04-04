from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import os
from Bio import PDB


files = ["5cep.pdb","5ceq.pdb"]


def target_sequence(file):
    i = 0
    with open(file) as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            if i > 0:
                raise TypeError("Dont put more than 1 sequence lol.")
            target = Seq(record.seq)
            i += 1
    return target


def psi_blast_results(file):
    list_sequences = []
    with open(file) as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            list_sequences.append(record)
    return list_sequences


def clustal_align(file):
    os.system(f"./clustalobin -i {file} -o output.sto --outfmt=st")
    return None


def load_structures(file_list):
    structures = []
    count = 1
    for file in file_list:
        Pdb_parser = PDB.PDBParser()
        struct_obj = Pdb_parser.get_structure(f"Structure{count}",file)
        structures.append(struct_obj)
    return structures

def superimpose(structure_list):
    x = PDB.Superimposer()
    x .set_atoms()
    pass


x = load_structures(files)
for structure in x:
    for atom in structure.get_atoms():
        print(atom)
