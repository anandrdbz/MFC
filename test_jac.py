import numpy as np 
import matplotlib.pyplot as plt 


#PERIODIC IN X, NEUMANN IN Y
N = 201
sol = np.zeros((N, N))
sol_old= np.zeros((N, N))
rhs = np.zeros((N, N)) + 0.01 * (1 + np.random.rand(N, N))
dx = (1 / N)
alpha = 10 * dx**2
rho = 10
N_its = 150
resid = np.zeros((N_its,))
omega = 1.6
for k in range(N_its):

	for j in range(N):
		for i in range(N):
			diag = (rho / alpha) + (4 / (dx**2))*rho
			if(j == 0):
				diag = diag - (1 / (dx**2))*rho
			if(j == 201):
				diag = diag - (1 / (dx**2))*rho
			sol[i, j] = rhs[i, j]
			if(j > 0):
				sol[i, j] = sol[i, j] + (1 / (dx**2))*rho*sol[i, j-1]
			if(j < N-1):
				sol[i, j] = sol[i, j] + (1 / (dx**2))*rho*sol[i, j+1]

			if(i == 0):
				sol[i, j] = sol[i,j] + (1 / (dx**2))*rho*sol[N-1, j]
			else: 
				sol[i, j] = sol[i, j] + (1 / (dx**2))*rho*sol[i-1, j]

			if(i == N-1):
				sol[i, j] = sol[i, j] + (1 / (dx**2))*rho*sol[0, j]
			else: 
				sol[i, j] = sol[i, j] + (1 / (dx**2))*rho*sol[i+1, j]

			sol[i, j]  = omega * (1/diag)*sol[i, j] + (1 - omega)*sol_old[i, j]

	
	resid[k] = np.max(sol - sol_old) / np.max(sol)
	print ("RESID", resid[k])

	for j in range(N):
		for i in range(N):
			sol_old[i, j] = 0 + sol[i, j]

plt.semilogy(np.arange(N_its) + 1, resid)
plt.xlabel("Iterations")
plt.ylabel("Residual")
plt.title("ALPHA =" + str(alpha / dx**2) + " dx^2")
plt.savefig("alpha_" + str(alpha / dx**2)  + ".png", dpi = 100)
plt.show()


