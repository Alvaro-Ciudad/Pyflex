import prody
import matplotlib.pyplot as plt


def write_pdb_CA(target, b_factor_array, type, path, name):
    """ Write CA-alpha Simulated B factors in pdb format"""
    protein = prody.parsePDB(f"{path}/temp/{target}.pdb")
    calphas = protein.select('calpha')
    calphas.setBetas(b_factor_array)
    prody.writePDB(f"{path}/outputs/{name}_{type}.pdb", calphas)


def plot_b_factors(b_factor_array, path, out, colour='orange'):
    """Plot B factors of CA against residue number"""
    plt.plot(b_factor_array, f'{colour}', label='Simulated')
    plt.grid()
    plt.legend()
    plt.savefig(f"{path}/outputs/{out}_scores_graph.png")
    plt.clf()
