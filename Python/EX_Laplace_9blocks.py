import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

from trAOSM import AOSM
from MAT_SymLaplace import SymLaplace
from DOM_Laplace import Laplace, calcError

def main():
    n = 99
    N = n ** 2

    # numerical parameters
    tol = 1e-8

    # domain
    block, trace, topRight, bottomLeft, rhsBlock, rhsTrace, uBlock, uTrace, indices, trace_indices = Laplace(N)

    solver = AOSM(block, trace, topRight, bottomLeft, rhsBlock, rhsTrace) # needs rewrite
    vBlock, vTrace = solver.main(tol, len(rhsTrace))

    # compare against the exact sparse solve
    err = calcError(uBlock, uTrace, vBlock, vTrace, indices, trace_indices, N)

if __name__ == "__main__":
    main()