import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

from trAOSM import AOSM
from MAT_SymLaplace import SymLaplace

def build_2blocks(n,h):
    x = np.linspace(-1, 1, 101)[1:-1]   # 99 interior points
    xx = np.tile(x, n)
    yy = np.repeat(x, n)

    ind = np.arange(n * n)

    cond1 = (xx < -h/2)
    cond2 = (xx > h/2)
    cond_trace = (xx>-h/2) & (xx<h/2)

    sub = [
        ind[cond1],
        ind[cond2]
    ]
    trace = ind[cond_trace]

    return sub, trace

def build_9blocks(n,h):
    x = np.linspace(-1, 1, 101)[1:-1]   # 99 interior points
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
    A, h = SymLaplace(N)
    f = -np.ones(N)

    # sub, trace = build_9blocks(n, h)
    ind_blocks, ind_trace = build_2blocks(n,h)

    blocks = [A[np.ix_(idx, idx)] for idx in ind_blocks]
    trace = A[np.ix_(ind_trace, ind_trace)]
    topRight = [A[np.ix_(idx, ind_trace)] for idx in ind_blocks]
    bottomLeft = [A[np.ix_(ind_trace, idx)] for idx in ind_blocks]
    rhsBlocks = [f[idx] for idx in ind_blocks]
    rhsTrace = f[ind_trace]

    solver = AOSM(blocks, trace, topRight, bottomLeft, rhsBlocks, rhsTrace)
    uBlocks, uTrace = solver.main()

    # assemble global solution roughly as in the MATLAB script
    u = np.zeros(N)
    for idx, ub in zip(ind_blocks, uBlocks):
        u[idx] = ub
    u[ind_trace] = uTrace

    # compare against the exact sparse solve
    u_exact = np.linalg.solve(A, f)
    err = np.linalg.norm(u - u_exact)

    print(f"Grid size: {n}x{n}")
    print(f"Trace size: {len(ind_trace)}")
    print(f"Global error norm: {err:.3e}")

        # Visualize
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    grid = int(np.sqrt(N))

    im1 = axes[0].imshow(u_exact.reshape(grid, grid), cmap="viridis")
    axes[0].set_title("Exact solution")
    fig.colorbar(im1, ax=axes[0], shrink=0.9)

    im2 = axes[1].imshow((u - u_exact).reshape(grid, grid), cmap="coolwarm")
    axes[1].set_title("Error: AOSM - exact")
    fig.colorbar(im2, ax=axes[1], shrink=0.9)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()