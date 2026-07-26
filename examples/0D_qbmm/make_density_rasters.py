from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

base = Path(".")

for i in range(1, 5):
    data = np.genfromtxt(
        base / f"density{i}.csv",
        delimiter=",",
        names=True,
    )

    xs = np.unique(data["R"])
    ys = np.unique(data["Rdot"])
    density = data["density"].reshape(len(ys), len(xs))

    fig = plt.figure(figsize=(4, 4), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])

    ax.imshow(
        density,
        origin="lower",
        extent=(xs.min(), xs.max(), ys.min(), ys.max()),
        cmap="Greys",
        norm=LogNorm(vmin=1.0e-4, vmax=1.0),
        interpolation="bilinear",
        aspect="auto",
    )

    ax.set_axis_off()

    fig.savefig(
        base / f"density{i}_raster.png",
        dpi=300,
        bbox_inches=None,
        pad_inches=0,
    )

    plt.close(fig)
