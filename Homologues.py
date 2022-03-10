from Bio import PDB
import warnings
class Homologues:
    def __init__(self, homologues_files, number=10):
        parser = PDB.PDBParser()
        self.homologues_objects = []
        try:
            homologues_files=homologues_files[:number]
        except IndexError:
            if len(homologues_files) < 10:
                warnings.warn("Number of homologues \
                inputed is less than recommended, whole list will be used")
            elif number > 10:
                warnings.warn("Incorrect input, 10 homologues will be used")
                homologues_files=homologues_files[:9]
        for file in homologues_files:
            try:
                struct_obj = parser.get_structure(f"{file}",file)
                self.homologues_objects.append(struct_obj)
            except:
                pass

    def __len__(self):
        return len(self.homologues_objects)

    def __str__(self):
        return f"Homologue with {len(self)} sequence objects."

    def get_sequences(self):
        return self.homologues_objects

    def get_CA_B_values_means(self):
        beta_factors = []
        for sequence in self.get_sequences():
            for residue in sequence.get_residues():
                try:
                    atom = residue["CA"]
                    beta_factors.append(atom.get_bfactor())
                except KeyError:
                    pass
        return sum(beta_factors)/len(beta_factors)




files = ["5cep.pdb","5ceq.pdb"]
x = Homologues(files, 3)
print(x.get_CA_B_values_means())
