import numpy as np

data = np.loadtxt("preston_final.csv", delimiter=",")
Nsave = 200

Rad_qbmm = data[:, 1]
Rad_mc = data[:, 4]

Pb_qbmm = data[:, 2]
Pb_mc = data[:, 5]

Mv_qbmm = data[:, 3]
Mv_mc = data[:, 6]

err_rad = (1 / np.sqrt(Nsave)) * np.sqrt(np.sum((1 - Rad_qbmm / Rad_mc) ** 2))
err_pb = (1 / np.sqrt(Nsave)) * np.sqrt(np.sum((1 - Pb_qbmm / Pb_mc) ** 2))
err_mv = (1 / np.sqrt(Nsave)) * np.sqrt(np.sum((1 - Mv_qbmm / Mv_mc) ** 2))

print(err_rad)
print(err_pb)
print(err_mv)
