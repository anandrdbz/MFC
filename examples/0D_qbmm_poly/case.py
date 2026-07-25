#!/usr/bin/env python3
import math
import json

# FLUID PROPERTIES
R_uni = 8314.0  # [J/kmol/K]
# Water
gam_l = 7.1  # [1]
pi_inf_l = 306.0e06  # [N/m2]
rho_l = 1.0e03  # [kg/m3]
mu_l = 1.002e-03  # [kg/m/s]
ss = 0.07275  # [kg/s2]
pv = 2.3388e03  # [N/m2]
vd = 0.242e-4  # [m2/s]

# Vapor
gam_v = 1.33  # [1]
M_v = 18.02  # [g/mol]
mu_v = 0.8816e-05  # [kg/m/s]
k_v = 0.019426  # [W/m/K]
cp_v = 2.1e3  # [J/kg/K]
R_v = R_uni / M_v  # [J/kg/K]

# Air
gam_g = 1.4  # [1]
M_g = 28.97  # [g/mol]
mu_g = 1.8e-05  # [kg/m/s]
k_g = 0.02556  # [W/m/K]
cp_g = 1.0e3  # [J/kg/K]
R_g = R_uni / M_g  # [J/kg/K]

# Bubble
R0ref = 10.0e-06  # [m]
p0ref = 101325.0  # [N/m2]
rho0ref = rho_l  # [kg/m3]
T0ref = 293.15  # [K]
vf0 = 0.00004  # [1]

# External condition
p0ext = 101325.0  # [N/m2]

# Pressure wave
pa = 1e5  # [N/m2]
freq = 300000.0  # [Hz]
cact = 1475.0  # [m/s]
wavelength = cact / freq  # [m]

# Reference values
x0 = R0ref  # [m]
p0 = p0ref  # [N/m2]
rho0 = rho0ref  # [kg/m3]
T0 = T0ref  # [K]
u0 = math.sqrt(p0 / rho0)  # [m/s]
t0 = x0 / u0  # [s]

#
cfl = 0.1
Nx = 400
Ldomain = 20.0e-03
L = Ldomain / x0
dt = 0.0005
Tfinal = 8
Nt = int(Tfinal / dt)
Nfiles = 200.0
Nout = int(math.ceil(Nt / Nfiles))
Nt = int(Nout * Nfiles)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -10.0e-03 / x0,
            "x_domain%end": 10.0e-03 / x0,
            "stretch_x": "F",
            "cyl_coord": "F",
            "m": 5,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": Nout,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 1,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "F",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "fd_order": 1,
            # 'schlieren_wrt'                :'T',
            "probe_wrt": "F",
            "num_probes": 1,
            "probe(1)%x": 0.0,
            # Patch 1 _ Background
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%length_x": 20.0e-03 / x0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": p0ext / p0,
            "patch_icpp(1)%alpha_rho(1)": (1.0 - vf0) * 1.0e03 / rho0,
            "patch_icpp(1)%alpha(1)": vf0,
            "patch_icpp(1)%r0": 1.0e000,
            "patch_icpp(1)%v0": 0.0e000,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_l - 1.0e00),
            "fluid_pp(1)%pi_inf": gam_l * (pi_inf_l / p0) / (gam_l - 1.0),
            # Bubbles
            "bubbles_euler": "T",
            "bub_0d": "T",
            "bubble_model": 3,
            "polytropic": "T",
            "polydisperse": "F",
            "thermal": 3,
            "nb": 1,
            # Bubble parameters
            "bub_pp%R0ref": R0ref / x0,
            "bub_pp%p0ref": p0ref / p0,
            "bub_pp%rho0ref": rho0ref / rho0,
            "bub_pp%ss": ss / (rho0 * x0 * u0 * u0),
            "bub_pp%pv": pv / p0,
            "bub_pp%mu_l": mu_l / (rho0 * x0 * u0),
            "bub_pp%gam_g": gam_g,
            "bub_pp%gam_v": gam_v,
            "bub_pp%T0ref": T0ref / T0,
            "bub_pp%vd": vd / (x0 * u0),
            "bub_pp%mu_v": mu_v / (rho0 * x0 * u0),
            "bub_pp%mu_g": mu_g / (rho0 * x0 * u0),
            "bub_pp%M_v": M_v,
            "bub_pp%M_g": M_g,
            "bub_pp%k_v": k_v * (T0 / (x0 * rho0 * u0 * u0 * u0)),
            "bub_pp%k_g": k_g * (T0 / (x0 * rho0 * u0 * u0 * u0)),
            "bub_pp%cp_v": cp_v * (T0 / (u0 * u0)),
            "bub_pp%cp_g": cp_g * (T0 / (u0 * u0)),
            "bub_pp%R_v": R_v * (T0 / (u0 * u0)),
            "bub_pp%R_g": R_g * (T0 / (u0 * u0)),
            "poly_sigma": 0.3,
            "qbmm": "T",
            "dist_type": 1,
            "sigR": 0.1,
            "sigV": 0.1,
            "rhoRV": 0.0,
        }
    )
)
