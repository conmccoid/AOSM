import numpy as np
import scipy.sparse as sp

def SymLaplace(N):
    # construct symmetric Laplacian matrix (2D) of size NxN for domain [-1,1]^2

    n = int(np.sqrt(N)) # NN = (n-1)^2
    h = 2/(n+1) # grid spacing
    d = np.ones(n)/h**2 # diagonal entries
    d = [d, -2*d, d]
    d = sp.diags(d, [-1, 0, 1], shape=(n, n)).toarray() # 1D Laplacian

    I = np.eye(n) # identity matrix
    A = np.kron(I, d) + np.kron(d, I) # 2D Laplacian
    return A, h

def sparseSymLaplace(N):
    # construct sparse symmetric Laplacian matrix (2D) of size NxN for domain [-1,1]^2

    n = int(np.sqrt(N)) # NN = (n-1)^2
    h = 2/(n+1) # grid spacing
    d = np.ones(n)/h**2 # diagonal entries
    d = [d, -2*d, d]
    d = sp.diags(d, [-1, 0, 1], shape=(n, n)) # 1D Laplacian

    I = sp.eye(n) # identity matrix
    A = sp.kron(I, d) + sp.kron(d, I) # 2D Laplacian
    return A