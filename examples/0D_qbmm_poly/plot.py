import numpy as np 
import matplotlib.pyplot as plt 

dt = 0.0005 
delt = 80
tstop = 16000
tstart = 0
count = 0
Nsamp = 1000

T_qbmm = []
Rad_qbmm = []
Rad_qbmm_pseudo_poly = []

Rad_moc = []

for i in range(tstart, tstop, delt):
	data_qbmm = np.loadtxt("D/prim.6.00." + str(i).zfill(6) + ".dat")[1,1]
	Rad_qbmm.append(data_qbmm)
	T_qbmm.append(count*delt*dt)

	data_qbmm = np.loadtxt("../0D_qbmm_pseudo_poly/D/prim.6.00." + str(i).zfill(6) + ".dat")[1,1]
	Rad_qbmm_pseudo_poly.append(data_qbmm)

	Rad_moc_temp = 0.0

	for j in range(1,Nsamp+1,1):
		data_moc = np.loadtxt("../0D_moc_poly/0D_moc_" + str(j) + "/D/prim.5.00." + str(i).zfill(6) + ".dat")[1,1]
		Rad_moc_temp = Rad_moc_temp + data_moc
	Rad_moc.append(Rad_moc_temp/Nsamp)

	count = count + 1
	print(i)

Nstep = np.size(Rad_qbmm)
poly_final = np.zeros((Nstep,4))
poly_final[:,0] = T_qbmm 
poly_final[:,1] = Rad_qbmm
poly_final[:,2] = Rad_moc
poly_final[:,3] = Rad_qbmm_pseudo_poly

np.savetxt("poly_final.csv", poly_final, delimiter = ",")

plt.plot(T_qbmm, Rad_qbmm, label = "qbmm")
plt.plot(T_qbmm, Rad_qbmm_pseudo_poly, label = "pseudo")
plt.plot(T_qbmm, Rad_moc, label = "Monte Carlo")
plt.savefig("../Rad_poly.png")
plt.close()


