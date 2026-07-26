import os
import random as rnd

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree

rnd.seed(1)


def fun_createDir(folderDir):
    isdir = os.path.isdir(folderDir)
    if isdir == False:
        os.mkdir(folderDir)


def rnd_log_normal(Rmean, Rdev):
    phi = np.sqrt(np.log(1.0 + (Rdev / Rmean) ** 2.0))
    mu = np.log(Rmean**2.0 / np.sqrt(Rmean**2.0 + Rdev**2.0))
    newRad = rnd.lognormvariate(mu, phi)
    return newRad


def cloud_generator_rectangular_prism(x_beg, x_end, y_beg, y_end, z_beg, z_end, voidFraction, Rmean, Rmin, Rmax, Rdev, Rskew, folderDir, Ncloud):

    volGas_target = voidFraction * (x_end - x_beg) * (y_end - y_beg) * (z_end - z_beg)
    volGas = 0.0
    Nbub_now = 0
    zeroVal = 0.0

    fileName = folderDir + "lag_bubbles.dat"
    file = open(fileName, "w")

    while True:
        # generate Radius
        newR = rnd_log_normal(Rmean, Rdev)

        # generate location
        new_x = rnd.uniform(x_beg, x_end)
        new_y = rnd.uniform(y_beg, y_end)
        new_z = rnd.uniform(z_beg, z_end)

        if newR >= Rmin and newR <= Rmax:
            Nbub_now = Nbub_now + 1
            volGas = volGas + (4.0 / 3.0) * np.pi * (newR**3.0)

            file.write("\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n" % (new_x, new_y, new_z, zeroVal, zeroVal, zeroVal, newR, zeroVal))

        if volGas >= volGas_target:
            break

    file.close()
    print(str(Ncloud) + " >> Cloud generated. Nbubs = " + str(Nbub_now) + " (" + fileName + ")")
    return Nbub_now


def nearest_neighbor(dataNumpy):
    tree = KDTree(dataNumpy)
    dist, ind = tree.query(dataNumpy, k=2)  # k=2: self + nearest neighbor
    nearest_distances = dist[:, 1]  # skip self (first one is always 0)
    mean_bubDist = np.mean(nearest_distances)
    min_bubDist = min(nearest_distances)
    return mean_bubDist, min_bubDist


def cloud_stats_plots_rectangular_prism(totVol, Rmean, gridSize, folderDir):
    filePath = folderDir + "lag_bubbles.dat"
    data_org = pd.read_csv(filePath, sep=r"\s+", header=None)
    data = data_org
    x = data[0]
    y = data[1]
    z = data[2]
    bubRad = data[6]
    outPath = folderDir + "summary.txt"
    f = open(outPath, "w")
    f.write(" >> Cloud generated. Nbubs = " + str(len(bubRad)) + " (" + filePath + ")\n")
    # Calculate final void fraction
    bubVol = 0.0
    for i in range(len(bubRad)):
        bubVol = bubVol + (4 / 3) * np.pi * bubRad[i] ** 3
    voidFraction = bubVol / totVol
    f.write("--> Calculated void fraction is %.6e \n" % voidFraction)
    dataNumpy = np.stack([x, y, z], axis=1).astype(np.float32)
    mean_bubDist, min_bubDist = nearest_neighbor(dataNumpy)
    f.write("Avg interbubble distance is %.4f times the mean bubble radius\n" % (mean_bubDist / Rmean))
    f.write("Minimum interbubble distance is %.4f times the mean bubble radius\n" % (min_bubDist / Rmean))
    f.write("Model inequality --> Grid size must be smaller than %.4f times the mean bubble radius.\n" % (mean_bubDist / Rmean))
    f.write("Model inequality --> Grid is %.4f times the mean bubble radius.\n\n" % (gridSize / Rmean))
    if gridSize < mean_bubDist:
        f.write("Model inequality passed! :)")
    else:
        f.write("Check inputs: Model inequality failed!!!!!!!! ")
    f.close()
    if gridSize < mean_bubDist:
        f = open("inequality_passed.txt", "a")
        f.write("%s, Nbubs:%i, mean and min bub-bub distance: %.4f & %.4f (> grid: %.4f)\n" % (folderDir, int(len(bubRad)), mean_bubDist, min_bubDist, gridSize))
        f.close()
    else:
        f = open("../inequality_failed.txt", "a")
        f.write("%s, Nbubs:%i, mean and min bub-bub distance: %.4f & %.4f (> grid: %.4f)\n" % (folderDir, int(len(bubRad)), mean_bubDist, min_bubDist, gridSize))
        f.close()
    # Plots
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(x, y, z, c=bubRad, cmap="jet")
    plt.colorbar(sc, label="Bubble Radius")
    name = folderDir + "plot_scatterDistribution.png"
    plt.savefig(name, transparent=False)
    plt.close()
    plt.hist(bubRad, density=True, bins=20)
    plt.xlabel("Radius of the bubbles")
    plt.ylabel("Normalized frequency")
    name = folderDir + "plot_sizeDistribution.png"
    plt.savefig(name, transparent=False)
    plt.close()
