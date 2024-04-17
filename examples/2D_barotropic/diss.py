import numpy as np
import matplotlib.pyplot as plt
arr = np.loadtxt("En_tot_igr.dat")
arr2 = np.loadtxt("En_tot_weno.dat")
plt.semilogy(np.arange(np.size(arr)), arr, label = "igr + lw")
plt.semilogy(np.arange(np.size(arr)), arr2, label = "weno5 + rk")
plt.legend()
plt.savefig("diss.png", dpi=300)
