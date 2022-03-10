import prody
import matplotlib.pyplot as plt


def GNMBuild(target):
    protein = prody.parsePDB(target)
    gnm = prody.GNM()
    calphas = protein.select('calpha')
    gnm.buildKirchhoff(calphas)
    gnm.calcModes(None)
    return gnm

def obtain_MSFs(target):
    network = GNMBuild(target)
    msfs = prody.calcSqFlucts(network)
    msfs_norm = msfs/prody.mean(msfs)
    return msfs_norm

def B_factor_scaling(target, homologues_obj):
    b_factor_mean = homologues_obj.get_CA_B_values_means()
    msfs_normalized = obtain_MSFs(target)
    return b_factor_mean*msfs_normalized

def write_pdb_CA(target, b_factor_array):
    protein = prody.parsePDB(target)
    name = target.split(".")
    calphas = protein.select('calpha')
    calphas.setBetas(b_factor_array)
    prody.writePDB(f"{name[0]}_simulated.pdb", calphas)

def plot_b_factors(b_factor_array):
    print(b_factor_array)
    plt.figure(figsize=(9, 5), dpi=300)
    plt.plot(b_factor_array, 'orange', label='Simulated')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    with open("temp2.pdb") as f:
        file = f.readlines()
    for index, line in enumerate(file):
        file[index] = line.strip()
    plot_b_factors(file)
