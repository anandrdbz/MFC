import numpy as np 

data = np.loadtxt('poly_final.csv', delimiter = ',')
Nsave = 200

Rad_qbmm  = data[:, 1]
Rad_mc = data[:, 2]
Rad_qbmm_pseudo = data[:, 3]

err_rad = (1/np.sqrt(Nsave)) * np.sqrt(np.sum((1 - Rad_qbmm/Rad_mc)**2))
err_rad_pseudo = (1/np.sqrt(Nsave)) * np.sqrt(np.sum((1 - Rad_qbmm_pseudo/Rad_mc)**2))

print(err_rad)
print(err_rad_pseudo)
