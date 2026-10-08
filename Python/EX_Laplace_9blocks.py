import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

from trAOSM import AOSM
from MAT_SymLaplace import SymLaplace
from DOM_Laplace import Laplace, calcError

def build_9blocks(n,h):
    x = np.linspace(-1, 1, n+2)[1:-1]   # n interior points
    xx = np.tile(x, n)
    yy = np.repeat(x, n)

    ind = np.arange(n * n)

    cond1 = ((xx < -0.3 -h/2) & (yy < -0.3-h/2)) | (np.isclose(xx, -0.3) & np.isclose(yy, -0.3))
    cond2 = (xx > -0.3 +h/2) & (xx < 0.3 -h/2) & (yy < -0.3)
    cond3 = ((xx > 0.3 +h/2) & (yy < -0.3-h/2)) | (np.isclose(xx, 0.3) & np.isclose(yy, -0.3))
    cond4 = (xx < -0.3 -h/2) & (yy > -0.3 +h/2) & (yy < 0.3 -h/2)
    cond5 = (xx > -0.3 +h/2) & (xx < 0.3 -h/2) & (yy > -0.3 +h/2) & (yy < 0.3 -h/2)
    cond6 = (xx > 0.3 +h/2) & (yy > -0.3 +h/2) & (yy < 0.3 -h/2)
    cond7 = ((xx < -0.3 -h/2) & (yy > 0.3 +h/2)) | (np.isclose(xx, -0.3) & np.isclose(yy, 0.3))
    cond8 = (xx > -0.3 +h/2) & (xx < 0.3 -h/2) & (yy > 0.3 +h/2)
    cond9 = ((xx > 0.3 +h/2) & (yy > 0.3 +h/2)) | (np.isclose(xx, 0.3) & np.isclose(yy, 0.3))

    trace_mask = (np.isclose(xx, -0.3) | np.isclose(xx, 0.3)) ^ (np.isclose(yy, -0.3) | np.isclose(yy, 0.3))

    sub = [
        ind[cond1],
        ind[cond2],
        ind[cond3],
        ind[cond4],
        ind[cond5],
        ind[cond6],
        ind[cond7],
        ind[cond8],
        ind[cond9]
    ]
    trace = ind[trace_mask]

    return sub, trace

def main():
    n = 99
    N = n ** 2

    # numerical parameters
    tol = 1e-8

    # domain
    block, trace, topRight, bottomLeft, rhsBlock, rhsTrace, uBlock, uTrace = Laplace(N)

    solver = AOSM(block, trace, topRight, bottomLeft, rhsBlock, rhsTrace) # needs rewrite
    vBlock, vTrace = solver.main(tol, len(rhsTrace))

    # compare against the exact sparse solve
    err = calcError(uBlock, uTrace, vBlock, vTrace, N)

if __name__ == "__main__":
    main()