from Bio import PDB
import sys
import requests


class Homologues:
    """ Homologue list processor.

    Needs a lists of files present in the directory/paths to be processed and
    the number of homologues that want to be saved as an optional argument.

    Includes functions to return length, return sequence objects, and calculate
    the mean of the B-factors of the Alpha Carbons.

    """

    def __init__(self, homologues_files):
        """ Class Init.

        Mandatory argument:
        homologues_files = list to be processed into a homologues object
        """

        parser = PDB.PDBParser()
        self.homologues_objects = []
        for file in homologues_files:
            file += ".pdb"
            try:
                struct_obj = parser.get_structure(f"{file}.pdb", file)
                self.homologues_objects.append(struct_obj)
            except FileNotFoundError:
                sys.stderr.write(
                    f"Warning, file not found: {file}. The file will not be processed.\n")
                pass

    def __len__(self):
        """Return length of homologues object as number of sequences inside"""
        return len(self.homologues_objects)

    def __str__(self):
        """String representation of the object as a function of its length"""
        return f"Homologue with {len(self)} sequence objects."

    def get_sequences(self):
        """Returns the list of sequence objects inside"""
        return self.homologues_objects

    def get_CA_B_values_means(self):
        """ B factor mean calculation.

        Computes the mean of the B factor values of the alpha carbons,
        of all the sequence objects inside.

         """
        beta_factors = []
        for sequence in self.get_sequences():
            for residue in sequence.get_residues():
                try:
                    atom = residue["CA"]
                    beta_factors.append(atom.get_bfactor())
                except KeyError:
                    pass
        return sum(beta_factors)/len(beta_factors)


def list_of_pdbfiles_preprocessor(file):
    file_list = []
    with open(file) as f:
        files = f.readlines()
    try:
        files = files[:10]
    except IndexError:
        sys.stderr.write(
            "Less than 10 homologues found, results might be worse.")
    for line in files:
        line = line.split(":")
        pdb = line[1]
        tup = pdb.split("_")
        tup[1] = tup[1].rstrip()
        file_list.append(tup)
    return file_list


def pdb_download_files(list_of_tuples_of_files):
    set_of_files = set(x[0] for x in list_of_tuples_of_files)
    for file in set_of_files:
        url = f"https://files.rcsb.org/download/{file}.pdb"
        r = requests.get(url, allow_redirects=True).content.decode("utf-8")
        with open(f"{file}.pdb", "w") as f:
            f.write(r)
    return None


def split_files_by_chain(list_of_tuples_of_files):
    list_of_homologues = ["_".join(x) for x in list_of_tuples_of_files]
    parser = PDB.PDBParser()
    io = PDB.PDBIO()
    set_of_files = set(x[0] for x in list_of_tuples_of_files)
    for file in set_of_files:
        structure = parser.get_structure(file, f"{file}.pdb")
        for chain in structure.get_chains():
            io.set_structure(chain)
            io.save(f"{structure.id}_{chain.id}.pdb")
    return list_of_homologues


def best_target_candidate(file_with_candidates):
    with open(file_with_candidates) as f:
        file = f.readlines()
    candidate1 = file[0].split(":")
    return candidate1[1]


if __name__ == '__main__':
    files = ["5cep.pdb", "5ceq.pdb", "lala.pdb"]
    #x = Homologues(files)
    #print(x.get_CA_B_values_means())
    #z = list_of_pdbfiles_preprocessor("files.txt")
    #print(pdb_download_files(z))
    print(best_target_candidate("best_hit.txt"))
