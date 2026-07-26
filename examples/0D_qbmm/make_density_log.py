#!/usr/bin/env python3
"""Prepare the PGFPlots input for Fig. 3 from the raw Monte Carlo density.

Reads the committed density{1..4}.csv -- the 220 x 220 regular grid over
(R, Rdot) written by convert_density.py, with values already clamped at
1e-4 -- and writes density{k}_log.csv, which quad_nodes.tex plots directly
as a native `surf`.

Two files are written per time:

  * density{k}_log.csv       -- full 220 x 220, used under LuaLaTeX
  * density{k}_log_coarse.csv -- 2x2 block-averaged 110 x 110, used under
                                 pdfLaTeX, whose fixed main memory cannot
                                 hold four 220 x 220 surf plots

The only other transformation is base-10 log of the density.  Doing it
here rather than as a `meta expr` in the .tex file is purely a
compile-time matter: evaluating ln() for 193,600 rows in TeX roughly
doubles the build (185 s vs 82 s with LuaLaTeX).  The grid, the ordering,
and the values are otherwise untouched, so density{k}.csv remains the
single source of truth.

The colorbar in quad_nodes.tex is labelled with the corresponding decades,
1e-4 to 1e0, matching `point meta min=-4` / `point meta max=0`.

Run from this directory:  python3 make_density_log.py
"""

import numpy as np

FLOOR = 1e-4  # matches the clamp applied when density{k}.csv was written
COARSE = 2  # 220 -> 110 for the pdfTeX fallback


def write_grid(path_out, R, V, logf):
    np.savetxt(
        path_out,
        np.column_stack([R.ravel(), V.ravel(), logf.ravel()]),
        delimiter=",",
        header="R,Rdot,logf",
        comments="",
        fmt="%.4f",
    )
    print(f"{path_out}: {logf.shape[0]}x{logf.shape[1]}, "
          f"log10 f in [{logf.min():.2f}, {logf.max():.2f}]")


def to_log(k):
    d = np.loadtxt(f"density{k}.csv", delimiter=",", skiprows=1)
    R = np.unique(d[:, 0])
    V = np.unique(d[:, 1])
    if len(R) * len(V) != d.shape[0]:
        raise ValueError(f"density{k}.csv: not a full {len(V)}x{len(R)} grid")

    Z = d[:, 2].reshape(len(V), len(R))  # rows written with R varying fastest
    RR, VV = np.meshgrid(R, V)
    write_grid(f"density{k}_log.csv", RR, VV, np.log10(np.clip(Z, FLOOR, None)))

    if len(R) % COARSE or len(V) % COARSE:
        raise ValueError(f"density{k}.csv: grid not divisible by {COARSE}")
    Zc = Z.reshape(len(V) // COARSE, COARSE, len(R) // COARSE, COARSE).mean(axis=(1, 3))
    Rc = R.reshape(-1, COARSE).mean(1)
    Vc = V.reshape(-1, COARSE).mean(1)
    RRc, VVc = np.meshgrid(Rc, Vc)
    write_grid(f"density{k}_log_coarse.csv", RRc, VVc,
               np.log10(np.clip(Zc, FLOOR, None)))


if __name__ == "__main__":
    for k in (1, 2, 3, 4):
        to_log(k)
