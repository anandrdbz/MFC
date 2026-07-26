import matplotlib.pyplot as plt
import numpy as np

N_end_qbmm = 1799
N_end_moc = 1799
N_end_mc = 1799

EL_results = np.zeros((N_end_qbmm, 4))
EL_extra_results = np.zeros((N_end_qbmm, 42))

data_qbmm = np.loadtxt("3D_bubblescreen_qbmm/D/probe1_prim.dat", usecols=(0, 5))
T_qbmm = data_qbmm[:N_end_qbmm, 0]
Pres_qbmm = data_qbmm[:N_end_qbmm, 1]
plt.plot(T_qbmm, Pres_qbmm, label="QBMM")
EL_results[:, 0] = T_qbmm
EL_results[:, 1] = Pres_qbmm.copy()

data_moc = np.loadtxt("3D_bubblescreen_moc/D/probe1_prim.dat", usecols=(0, 5))
T_moc = data_moc[:N_end_moc, 0]
Pres_moc = data_moc[:N_end_moc, 1]
plt.plot(T_moc, Pres_moc, label="Method of Classes")
EL_results[:, 2] = Pres_moc.copy()

count = 0
for i in range(1, 41):
    data_mc = np.loadtxt("3D_lagrange_qbmm_" + str(i) + "/D/probe1_prim.dat", usecols=(0, 5))
    if count == 0:
        T_mc = data_mc[:N_end_mc, 0]
        Pres_mc = 0.0 * T_mc
    Pres_mc = Pres_mc + data_mc[:N_end_mc, 1]
    EL_extra_results[:, count] = data_mc[:N_end_mc, 1].copy()
    count = count + 1
    print(i)

Pres_mc = Pres_mc / count
EL_results[:, 3] = Pres_mc.copy()


l2 = np.sqrt(np.sum((1 - Pres_qbmm / Pres_mc) ** 2)) / np.sqrt(N_end_qbmm)
print(l2)

# np.savetxt("EL_results.csv", EL_results, delimiter = ",")
# np.savetxt("EL_extra_results.csv", EL_extra_results, delimiter = ",")

# plt.plot(T_mc, Pres_mc, label = "Euler--Lagrange (Mean)")
# plt.legend()
# plt.xlabel("t/t0")
# plt.ylabel("p/p0")
# plt.xlim(2, 25)
# plt.savefig("final.png", dpi = 200)
