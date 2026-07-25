from py_functions import *
import math

# Reference scaling
x0 = 10e-6  # reference length (m)

# Domain (normalized by x0)
H = 5.e-3/ x0  # cubic cloud size

x_beg, x_end = -0.5 * H, 0.5 * H 
y_beg, y_end = -0.5 * H, 0.5 * H
z_beg, z_end = -0.5 * H, 0.5 * H

# Bubble / cloud properties
void_fraction = 4e-5

Rmean = 10e-6 / x0  # mean radius

# Lognormal dispersion (monodisperse: Rdev_log = 0.0)
Rdev_log = 0.2
Rdev = Rmean * math.sqrt(math.exp(Rdev_log**2) - 1.0)

# Min/max bounds (from MFC common/m_helper.fpp -> s_simpson)
Rmin = Rmean * (0.8 * math.exp(-2.8 * Rdev_log))
Rmax = Rmean * (0.2 * math.exp(9.5 * Rdev_log) + 1.0)

Rskew = 0.0

# Grid resolution
grid_size = 100e-6 / x0

# Output
output_path = "input/"
fun_createDir(output_path)

# Generate cloud
cloud_volume = (x_end - x_beg) * (y_end - y_beg) * (z_end - z_beg)

n_bubbles = cloud_generator_rectangular_prism(
                x_beg, x_end, y_beg, y_end, z_beg, z_end,
                void_fraction, Rmean, Rmin, Rmax, Rdev_log, Rskew,
                output_path, 1)

cloud_stats_plots_rectangular_prism( cloud_volume, Rmean, grid_size, output_path)
