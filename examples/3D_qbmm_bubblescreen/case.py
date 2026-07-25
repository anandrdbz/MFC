import json
import math
import os

# Reference values for nondimensionalization
x0 = 10e-6
rho0 = 1e3
p0 = 101325
c0 = math.sqrt(p0 / rho0)
T0 = 298

r0 = 10e-6 # Mean bubble radius

# Domain bounds and resolution (w/ stretching)
num_cells = [400, 50, 50]
dim = [x * r0 / x0 for x in [4000.0, 500.0, 500.0]]

# Water properties
rho_w = 1000.0 / rho0
gamma_w = 6.12
pi_inf_w = 3.43e8 / p0
mu_w = 1e-3 / (rho0 * c0 * x0)

# Gas properties
gamma_g = 1.4
mu_g = 1.48e-5 / (rho0 * c0 * x0)

# Lagrangian bubble properties
R_uni = 8314          # Universal gas constant - J/kmol/K
MW_g = 28.0           # Molar weight of the gas - kg/kmol
MW_v = 18.0           # Molar weight of the vapor - kg/kmol
gam_g = 1.4           # Specific heat ratio of the gas
gam_v = 1.333         # Specific heat ratio of the vapor
pv = 2350             # Vapor pressure of the host - Pa
cp_g = 1.0e3          # Specific heat of the gas - J/kg/K
cp_v = 2.1e3          # Specific heat of the vapor - J/kg/K
k_g = 0.025           # Thermal conductivity of the gas - W/m/K
k_v = 0.02            # Thermal conductivity of the vapor - W/m/K
diffVapor = 2.5e-5    # Diffusivity coefficient of the vapor - m2/s
sigBubble = 0.074     # Surface tension of the bubble - N/m

# Ambient properties
pAmb = 101325 / p0

# Sound speed and delta t
dx = dim[0] / num_cells[0]
cw = math.sqrt(gamma_w * (pAmb + pi_inf_w) / rho_w)

# Acoustic source params
pAc = 1e5 / p0
pulse_width_factor = 2  # Multiplier on Gaussian pulse width

# Minnaert (natural) frequency for mean bubble radius (nondim R0=1)
f_nat = (1.0 / (2 * math.pi)) * math.sqrt(3 * gamma_g * pAmb / rho_w)
lambda_nat = cw / f_nat

f_ac = 0.3 
lambda_ac = cw / f_ac

# Timestepping and output parameters
r0 = r0 / x0  # Nondimensional reference radius
T_RC = 0.915 * r0 * math.sqrt(rho_w / (pAmb - pv / p0))
n_collapse_times = 27
nsave = 100
tend = n_collapse_times * T_RC

tsave = tend / nsave
cfl = 0.2
eps = 1.e-7

dt = cfl * dx / cw
t_step_start = 0
t_step_stop = int(tend / dt) + 1
t_step_save = int(t_step_stop / nsave)

vf0 = 4.e-5

case = {
    "run_time_info": "T",
    "parallel_io": "T",
    "probe_wrt": "T",
    "fd_order": 2,
    "num_probes": 1,
    "probe(1)%x": 0.0,
    "probe(1)%y": 0.0,
    "probe(1)%z": 0.0,
    "rdma_mpi": "F",
    "m": num_cells[0]-1,
    "n": num_cells[1]-1,
    "p": num_cells[2]-1,
    "dt": dt,
    "t_step_start": t_step_start,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "num_patches": 2,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "recon_type": 1,
    "weno_order": 3,
    "mapped_weno": "T",
    "weno_eps": 1e-16,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "viscous": "T",
    "precision": 2,
    "format": 1,
    "prim_vars_wrt": "T",
    "lag_db_wrt": "T",
    "lag_txt_wrt": "T",
    "x_domain%beg": -0.5*dim[0],
    "x_domain%end": 0.5*dim[0],
    "y_domain%beg": -0.5*dim[1],
    "y_domain%end": 0.5*dim[1],
    "z_domain%beg": -0.5*dim[2],
    "z_domain%end": 0.5*dim[2],
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    # Background flow: pure water
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_z": dim[2],
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_y": dim[1],
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%length_x": dim[0],
    "patch_icpp(1)%pres": pAmb,
    "patch_icpp(1)%alpha_rho(1)": (1 - eps)*rho_w,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%r0": 1.0,
    "patch_icpp(1)%v0": 0.0e00,
    #Bubble Screen 
    "patch_icpp(2)%geometry": 9,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%length_z": 500.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%length_y": 500.0,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%length_x": 500.0,
    "patch_icpp(2)%pres": pAmb,
    "patch_icpp(2)%alpha_rho(1)": (1 - vf0)*rho_w,
    "patch_icpp(2)%alpha(1)": vf0,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%r0": 1.0,
    "patch_icpp(2)%v0": 0.0e00,
    "patch_icpp(2)%alter_patch(1)": "T",
    # Fluid parameters
    # Fluid 1: Water (host medium)
    "fluid_pp(1)%gamma": 1.0 / (gamma_w - 1.0),
    "fluid_pp(1)%pi_inf": gamma_w * pi_inf_w / (gamma_w - 1.0),
    "fluid_pp(1)%Re(1)": 1.0 / mu_w,
    #Euler Bubbles
    "bubbles_euler": "T",
    "bubble_model": 2,              # Keller-Miksis model
    "thermal": 3,
    "polytropic": "F",
    "nb": 51,
    "poly_sigma": 0.3,
    "qbmm": "T",
    "dist_type": 1,
    "sigR": 0.2,
    "sigV": 0.2,
    "rhoRV": 0.0,
    "adap_dt": "F",
    # Bubble parameters (nondimensionalized)
    "bub_pp%R0ref": 1.0,
    "bub_pp%p0ref": 1.0,
    "bub_pp%rho0ref": 1.0,
    "bub_pp%T0ref": 1.0,
    "bub_pp%ss": sigBubble / (rho0 * x0 * c0 * c0),
    "bub_pp%pv": pv / p0,
    "bub_pp%vd": diffVapor / (x0 * c0),
    "bub_pp%mu_l": mu_w,
    "bub_pp%gam_v": gam_v,
    "bub_pp%gam_g": gam_g,
    "bub_pp%M_v": MW_v,
    "bub_pp%M_g": MW_g,
    "bub_pp%k_v": k_v * (T0 / (x0 * rho0 * c0 * c0 * c0)),
    "bub_pp%k_g": k_g * (T0 / (x0 * rho0 * c0 * c0 * c0)),
    "bub_pp%cp_v": cp_v * (T0 / (c0 * c0)),
    "bub_pp%cp_g": cp_g * (T0 / (c0 * c0)),
    "bub_pp%R_v": (R_uni / MW_v) * (T0 / (c0 * c0)),
    "bub_pp%R_g": (R_uni / MW_g) * (T0 / (c0 * c0)),
    # Acoustic source
    "acoustic_source": "T",
    "num_source": 1,
    "acoustic(1)%support": 3,
    "acoustic(1)%npulse": 1,
    "acoustic(1)%mag": pAc,
    "acoustic(1)%pulse": 1,
    "acoustic(1)%wavelength" : lambda_ac,
    "acoustic(1)%length": 2*dim[2],
    "acoustic(1)%height": 2*dim[1],
    "acoustic(1)%loc(1)": -700.0,
    "acoustic(1)%loc(2)": 0.0,
    "acoustic(1)%loc(3)": 0.0,
    "acoustic(1)%dir": 0.0,
    "acoustic(1)%delay": 0.0,
}

print(json.dumps(case, indent=4))
