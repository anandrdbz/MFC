import matplotlib.pyplot as plt
import numpy as np

SGM_EPS = 10**-16


def s_hyqmom(fmom):
    """
    Symmetric two-node quadrature based on moments [m0, m1, m2].
    """
    fmom = np.asarray(fmom, dtype=float)

    bu = fmom[1] / fmom[0]
    d2 = fmom[2] / fmom[0]
    c2 = max(d2 - bu**2, SGM_EPS)

    frho = np.array(
        [
            fmom[0] / 2.0,
            fmom[0] / 2.0,
        ]
    )

    fup = np.array(
        [
            bu - np.sqrt(c2),
            bu + np.sqrt(c2),
        ]
    )

    return frho, fup


def compute_2d_quadrature(momin):
    """
    Construct a four-node conditional quadrature from the six moments
    [m00, m10, m01, m20, m11, m02].
    """
    momin = np.asarray(momin, dtype=float)

    m00, m10, m01, m20, m11, m02 = momin

    bu = m10 / m00
    bv = m01 / m00

    d20 = m20 / m00
    d11 = m11 / m00
    d02 = m02 / m00

    c20 = d20 - bu**2
    c11 = d11 - bu * bv
    c02 = d02 - bv**2

    # X quadrature
    M1 = np.array([1.0, 0.0, c20])
    myrho, up = s_hyqmom(M1)

    # Conditional mean in Y
    Vf = c11 * up / c20

    # Conditional variance in Y
    mu2 = max(
        SGM_EPS,
        c02 - np.sum(myrho * Vf**2),
    )

    # Y quadrature
    M3 = np.array([1.0, 0.0, mu2])
    myrho3, up3 = s_hyqmom(M3)

    # Weights
    wght = m00 * np.array(
        [
            myrho[0] * myrho3[0],
            myrho[0] * myrho3[1],
            myrho[1] * myrho3[0],
            myrho[1] * myrho3[1],
        ]
    )

    # Abscissas
    abscX = bu + np.array(
        [
            up[0],
            up[0],
            up[1],
            up[1],
        ]
    )

    abscY = bv + np.array(
        [
            Vf[0] + up3[0],
            Vf[0] + up3[1],
            Vf[1] + up3[0],
            Vf[1] + up3[1],
        ]
    )

    return wght, abscX, abscY


dt = 0.0005
delt = 80
tstop = 16000
tstart = 0
count = 0
Nsamp = 1000

T_qbmm = []
Rad_qbmm = []
Pb_qbmm = []
Mv_qbmm = []

Rad_moc = []
Pb_moc = []
Mv_moc = []

Rad_mc = np.zeros((Nsamp, 4))
Rdot_mc = np.zeros((Nsamp, 4))

momin = np.zeros(6)
node_locs = np.zeros((4, 2))

for i in range(tstart, tstop, delt):
    data_qbmm = np.loadtxt("D/prim.6.00." + str(i).zfill(6) + ".dat")[1, 1]
    Rad_qbmm.append(data_qbmm)
    T_qbmm.append(count * delt * dt)

    if count == 32:
        momin[0] = 1.0
        momin[1] = data_qbmm
        momin[2] = np.loadtxt("D/prim.7.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[3] = np.loadtxt("D/prim.8.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[4] = np.loadtxt("D/prim.9.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[5] = np.loadtxt("D/prim.10.00." + str(i).zfill(6) + ".dat")[1, 1]

        wght, abscX, abscY = compute_2d_quadrature(momin)
        node_locs[:, 0] = abscX
        node_locs[:, 1] = abscY

        np.savetxt("Node1.csv", node_locs, delimiter=",")

    if count == 62:
        momin[0] = 1.0
        momin[1] = data_qbmm
        momin[2] = np.loadtxt("D/prim.7.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[3] = np.loadtxt("D/prim.8.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[4] = np.loadtxt("D/prim.9.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[5] = np.loadtxt("D/prim.10.00." + str(i).zfill(6) + ".dat")[1, 1]

        wght, abscX, abscY = compute_2d_quadrature(momin)
        node_locs[:, 0] = abscX
        node_locs[:, 1] = abscY

        np.savetxt("Node2.csv", node_locs, delimiter=",")

    if count == 92:
        momin[0] = 1.0
        momin[1] = data_qbmm
        momin[2] = np.loadtxt("D/prim.7.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[3] = np.loadtxt("D/prim.8.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[4] = np.loadtxt("D/prim.9.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[5] = np.loadtxt("D/prim.10.00." + str(i).zfill(6) + ".dat")[1, 1]

        wght, abscX, abscY = compute_2d_quadrature(momin)
        node_locs[:, 0] = abscX
        node_locs[:, 1] = abscY

        np.savetxt("Node3.csv", node_locs, delimiter=",")

    if count == 162:
        momin[0] = 1.0
        momin[1] = data_qbmm
        momin[2] = np.loadtxt("D/prim.7.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[3] = np.loadtxt("D/prim.8.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[4] = np.loadtxt("D/prim.9.00." + str(i).zfill(6) + ".dat")[1, 1]
        momin[5] = np.loadtxt("D/prim.10.00." + str(i).zfill(6) + ".dat")[1, 1]

        wght, abscX, abscY = compute_2d_quadrature(momin)
        node_locs[:, 0] = abscX
        node_locs[:, 1] = abscY

        np.savetxt("Node4.csv", node_locs, delimiter=",")

    data_qbmm = 0.0
    for q in range(1, 5, 1):
        data_qbmm = data_qbmm + np.loadtxt("D/pres.1." + str(q) + ".00." + str(i).zfill(6) + ".dat")[1, 1]
    Pb_qbmm.append(data_qbmm / 4.0)

    data_qbmm = 0.0
    for q in range(1, 5, 1):
        data_qbmm = data_qbmm + np.loadtxt("D/mv.1." + str(q) + ".00." + str(i).zfill(6) + ".dat")[1, 1]
    Mv_qbmm.append(data_qbmm / 4.0)

    Rad_moc_temp = 0.0
    Pb_moc_temp = 0.0
    Mv_moc_temp = 0.0

    for j in range(1, Nsamp + 1, 1):
        data_moc = np.loadtxt("../0D_moc/0D_moc_" + str(j) + "/D/prim.5.00." + str(i).zfill(6) + ".dat")[1, 1]
        Rad_moc_temp = Rad_moc_temp + data_moc

        if count == 32:
            Rad_mc[j - 1, 0] = data_moc
        elif count == 62:
            Rad_mc[j - 1, 1] = data_moc
        elif count == 92:
            Rad_mc[j - 1, 2] = data_moc
        elif count == 162:
            Rad_mc[j - 1, 3] = data_moc

        data_moc = np.loadtxt("../0D_moc/0D_moc_" + str(j) + "/D/prim.6.00." + str(i).zfill(6) + ".dat")[1, 1]

        if count == 32:
            Rdot_mc[j - 1, 0] = data_moc
        elif count == 62:
            Rdot_mc[j - 1, 1] = data_moc
        elif count == 92:
            Rdot_mc[j - 1, 2] = data_moc
        elif count == 162:
            Rdot_mc[j - 1, 3] = data_moc

    Rad_moc.append(Rad_moc_temp / Nsamp)

    for j in range(1, Nsamp + 1, 1):
        data_moc = np.loadtxt("../0D_moc/0D_moc_" + str(j) + "/D/prim.7.00." + str(i).zfill(6) + ".dat")[1, 1]
        Pb_moc_temp = Pb_moc_temp + data_moc
    Pb_moc.append(Pb_moc_temp / Nsamp)

    for j in range(1, Nsamp + 1, 1):
        data_moc = np.loadtxt("../0D_moc/0D_moc_" + str(j) + "/D/prim.8.00." + str(i).zfill(6) + ".dat")[1, 1]
        Mv_moc_temp = Mv_moc_temp + data_moc
    Mv_moc.append(Mv_moc_temp / Nsamp)

    count = count + 1
    print(i)

Nstep = np.size(Rad_qbmm)
preston_final = np.zeros((Nstep, 7))
preston_final[:, 0] = T_qbmm
preston_final[:, 1] = Rad_qbmm
preston_final[:, 2] = Pb_qbmm
preston_final[:, 3] = Mv_qbmm
preston_final[:, 4] = Rad_moc
preston_final[:, 5] = Pb_moc
preston_final[:, 6] = Mv_moc

np.savetxt("preston_final.csv", preston_final, delimiter=",")
np.savetxt("Rad_mc.csv", Rad_mc, delimiter=",")
np.savetxt("Rdot_mc.csv", Rdot_mc, delimiter=",")


plt.plot(T_qbmm, Rad_qbmm, label="qbmm")
plt.plot(T_qbmm, Rad_moc, label="Monte Carlo")
plt.savefig("../Rad.png")
plt.close()

plt.plot(T_qbmm, Pb_qbmm, label="qbmm")
plt.plot(T_qbmm, Pb_moc, label="Monte Carlo")
plt.savefig("../Pb.png")
plt.close()
