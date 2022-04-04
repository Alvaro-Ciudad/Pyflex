import prody
import numpy as np


def GNMBuild(target):
    """Obtain GNM object with modes calculated for alpha carbons from target"""
    protein = prody.parsePDB(target)
    gnm = prody.GNM()
    calphas = protein.select('calpha')
    gnm.buildKirchhoff(calphas)
    gnm.calcModes(None)
    return gnm


def obtain_MSFs(target):
    """Obtain Normalized Mean Square Fluctuations from target protein"""
    network = GNMBuild(target)
    msfs = prody.calcSqFlucts(network)
    msfs_norm = msfs/np.mean(msfs)
    return msfs_norm


def B_factor_scaling(homologues_obj, msfs_normalized):
    """Obtain simulated B factor from homologues and target"""
    b_factor_mean = homologues_obj.get_CA_B_values_means()
    return b_factor_mean*msfs_normalized


if __name__ == '__main__':
    with open("temp2.pdb") as f:
        file = f.readlines()
    for index, line in enumerate(file):
        file[index] = float(line.strip())
    #plot_b_factors(file)
