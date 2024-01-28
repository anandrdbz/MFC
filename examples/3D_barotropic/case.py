#!/usr/bin/env python3

import math
import json

# Numerical setup
Nx = 100
dx = 1./(1.*(Nx+1))

Tend = 0.2
maxwavespeed = 1.
lam = 0.2

mydt = lam*dx/maxwavespeed
Nt = int(Tend/mydt)
#Nt = 2
# print('T',Tend)

gam = 1.4

rho = 0.1
u = 1.

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'F',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : 0.E+00,
    'x_domain%end'                 : 1.E+00,
    'y_domain%beg'                 : 0.E+00,
    'y_domain%end'                 : 1.E+00,
    'z_domain%beg'                 : 0.E+00,
    'z_domain%end'                 : 1.E+00,
    'm'                            : Nx,
    'n'                            : Nx,
    'p'                            : Nx,
    'dt'                           : mydt,
    't_step_start'                 : 0,
    't_step_stop'                  : 3000,
    't_step_save'                  : 60,
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 2,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 4,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',  
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -1,
    'bc_x%end'                     : -1,
    'bc_y%beg'                     : -1,
    'bc_y%end'                     : -1,
    'bc_z%beg'                     : -1,
    'bc_z%end'                     : -1,
    'barotropic'                   : 'T',
    'cu_mpi'                   : 'F',
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================================


    'patch_icpp(1)%geometry'       : 9,
    'patch_icpp(1)%x_centroid'     : 0.25,
    'patch_icpp(1)%y_centroid'     : 0.5,
    'patch_icpp(1)%z_centroid'     : 0.5,
    'patch_icpp(1)%length_x'       : 0.5,
    'patch_icpp(1)%length_y'       : 1.0,
    'patch_icpp(1)%length_z'       : 1.0,    
    'patch_icpp(1)%vel(1)'         : 0.0*u,
    'patch_icpp(1)%vel(2)'         : 0.0*u,
    'patch_icpp(1)%vel(3)'         : 0.0*u,
    'patch_icpp(1)%alpha_rho(1)'   : rho, 
    'patch_icpp(1)%pres'           : 0.2*(rho)**1.4,
    'patch_icpp(1)%alpha(1)'       : 1.,

    'patch_icpp(2)%geometry'       : 9,
    'patch_icpp(2)%x_centroid'     : 0.75,
    'patch_icpp(2)%y_centroid'     : 0.5,
    'patch_icpp(2)%z_centroid'     : 0.5,
    'patch_icpp(2)%length_x'       : 0.5,
    'patch_icpp(2)%length_y'       : 1.0,
    'patch_icpp(2)%length_z'       : 1.0,
    'patch_icpp(2)%vel(1)'         : 0.0*u,
    'patch_icpp(2)%vel(2)'         : 0.0*u,
    'patch_icpp(2)%vel(3)'         : 0.0*u,
    'patch_icpp(2)%alpha_rho(1)'   : rho, 
    'patch_icpp(2)%pres'           : 0.2*(rho)**1.4,
    'patch_icpp(2)%alpha(1)'       : 1.,
                                                            

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.4,
    # 'fluid_pp(1)%gamma'            : 1.E+00/(gam-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.0,
    # ==========================================================================
}))

# ==============================================================================
