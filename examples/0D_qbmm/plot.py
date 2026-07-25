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
Pb_qbmm = []
Mv_qbmm = []

Rad_moc = []
Pb_moc = []
Mv_moc = []

Rad_mc = np.zeros((Nsamp, 4))
Rdot_mc = np.zeros((Nsamp, 4))

for i in range(tstart, tstop, delt):
	data_qbmm = np.loadtxt("D/prim.6.00." + str(i).zfill(6) + ".dat")[1,1]
	Rad_qbmm.append(data_qbmm)
	T_qbmm.append(count*delt*dt)

	data_qbmm = 0.0
	for q in range(1, 5, 1):
		data_qbmm = data_qbmm + np.loadtxt("D/pres.1."  + str(q) + ".00."+ str(i).zfill(6) + ".dat")[1,1]
	Pb_qbmm.append(data_qbmm/4.0)

	data_qbmm = 0.0
	for q in range(1, 5, 1):
		data_qbmm = data_qbmm + np.loadtxt("D/mv.1."  + str(q) + ".00."+ str(i).zfill(6) + ".dat")[1,1]
	Mv_qbmm.append(data_qbmm/4.0)


	Rad_moc_temp = 0.0
	Pb_moc_temp = 0.0
	Mv_moc_temp = 0.0

	for j in range(1,Nsamp+1,1):
		data_moc = np.loadtxt("../0D_moc/0D_moc_" + str(j) + "/D/prim.5.00." + str(i).zfill(6) + ".dat")[1,1]
		Rad_moc_temp = Rad_moc_temp + data_moc

		if(count == 32):
			Rad_mc[j-1,0] = data_moc 
		elif(count == 62):
			Rad_mc[j-1,1] = data_moc
		elif(count == 92):
			Rad_mc[j-1,2] = data_moc
		elif(count == 162):
			Rad_mc[j-1,3] = data_moc

		data_moc = np.loadtxt("../0D_moc/0D_moc_" + str(j) + "/D/prim.6.00." + str(i).zfill(6) + ".dat")[1,1]

		if(count == 32):
			Rdot_mc[j-1,0] = data_moc 
		elif(count == 62):
			Rdot_mc[j-1,1] = data_moc
		elif(count == 92):
			Rdot_mc[j-1,2] = data_moc
		elif(count == 162):
			Rdot_mc[j-1,3] = data_moc




	Rad_moc.append(Rad_moc_temp/Nsamp)

	for j in range(1,Nsamp+1,1):
		data_moc = np.loadtxt("../0D_moc/0D_moc_" + str(j) + "/D/prim.7.00." + str(i).zfill(6) + ".dat")[1,1]
		Pb_moc_temp = Pb_moc_temp + data_moc
	Pb_moc.append(Pb_moc_temp/Nsamp)

	for j in range(1,Nsamp+1,1):
		data_moc = np.loadtxt("../0D_moc/0D_moc_" + str(j) + "/D/prim.8.00." + str(i).zfill(6) + ".dat")[1,1]
		Mv_moc_temp = Mv_moc_temp + data_moc
	Mv_moc.append(Mv_moc_temp/Nsamp)

	count = count + 1
	print(i)

Nstep = np.size(Rad_qbmm)
preston_final = np.zeros((Nstep,7))
preston_final[:,0] = T_qbmm 
preston_final[:,1] = Rad_qbmm
preston_final[:,2] = Pb_qbmm
preston_final[:,3] = Mv_qbmm
preston_final[:,4] = Rad_moc
preston_final[:,5] = Pb_moc
preston_final[:,6] = Mv_moc

np.savetxt("preston_final.csv", preston_final, delimiter = ",")
np.savetxt("Rad_mc.csv", Rad_mc, delimiter = ",")
np.savetxt("Rdot_mc.csv", Rdot_mc, delimiter = ",")


plt.plot(T_qbmm, Rad_qbmm, label = "qbmm")
plt.plot(T_qbmm, Rad_moc, label = "Monte Carlo")
plt.savefig("../Rad.png")
plt.close()

plt.plot(T_qbmm, Pb_qbmm, label = "qbmm")
plt.plot(T_qbmm, Pb_moc, label = "Monte Carlo")
plt.savefig("../Pb.png")
plt.close()



