from pathlib import Path

import numpy as np
from scipy.stats import gaussian_kde

base = Path(".")

# Each file is assumed to contain 4 columns, one for each simulation time.
rad = np.loadtxt(base / "Rad_mc.csv", delimiter=",")
rdot = np.loadtxt(base / "Rdot_mc.csv", delimiter=",")

# Grid used for the PGFPlots contour data.
x = np.linspace(0.6, 1.4, 220)
y = np.linspace(-0.6, 0.6, 220)

X, Y = np.meshgrid(x, y)

query_points = np.vstack(
    [
        X.ravel(),
        Y.ravel(),
    ]
)

densities = []

for i in range(4):
    samples = np.vstack(
        [
            rad[:, i],
            rdot[:, i],
        ]
    )

    kde = gaussian_kde(samples)

    density = kde(query_points).reshape(X.shape)

    densities.append(density)

# Use one common normalization across all four panels.
density_max = max(density.max() for density in densities)

for i, density in enumerate(densities, start=1):
    density_normalized = density / density_max

    # Prevent zero values on the logarithmic color scale.
    density_normalized = np.maximum(
        density_normalized,
        1.0e-4,
    )

    output_file = base / f"density{i}.csv"

    with output_file.open("w") as file:
        file.write("R,Rdot,density\n")

        for iy in range(len(y)):
            for ix in range(len(x)):
                file.write(f"{X[iy, ix]:.12e},{Y[iy, ix]:.12e},{density_normalized[iy, ix]:.12e}\n")

            # Blank line separates adjacent grid rows for PGFPlots.
            file.write("\n")

    print(f"Created {output_file}")
