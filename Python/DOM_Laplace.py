import numpy as np
import scipy.sparse as sp

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nrank= comm.Get_size()

def Laplace(N):
    # dense Laplacian
    n = int(np.sqrt(N))
    h = 2/(n+1)
    if rank==0:
        d = np.ones(n)
        d = [d, -2*d, d]
        d = sp.diags(d, [-1,0,1], shape=(n,n)).toarray()
        I = np.eye(n)
        A = np.kron(I,d) + np.kron(d,I)
        f = -(h**2)*np.ones(N)
        u = np.linalg.solve(A,f)
        x = np.linspace(-1, 1, n+2)[1:-1] # n interior points
        xx = np.tile(x,n)
        # yy = np.repeat(x,n)
        ind = np.arange(N)
        if nrank==1:
            trace_mask = np.zeros(N, dtype=bool)
            block_mask = ~trace_mask
        elif nrank==2:
            trace_mask = (xx>-h/2) & (xx<h/2)
            block_mask = (xx < -h/2)
            other_mask = (xx > h/2)
            comm.send([
                A[np.ix_(ind[other_mask], ind[other_mask])],
                A[np.ix_(ind[other_mask], ind[trace_mask])],
                A[np.ix_(ind[trace_mask], ind[other_mask])],
                f[ind[other_mask]],
                u[ind[other_mask]],
            ], dest=1)
        trace = A[np.ix_(ind[trace_mask], ind[trace_mask])]
        rhsTrace = f[ind[trace_mask]]
        uTrace = u[ind[trace_mask]]
        trace_data = (trace, rhsTrace, uTrace)
        
        block = A[np.ix_(ind[block_mask],ind[block_mask])]
        topRight = A[np.ix_(ind[block_mask], ind[trace_mask])]
        bottomLeft = A[np.ix_(ind[trace_mask], ind[block_mask])]
        rhsBlock = f[ind[block_mask]]
        uBlock = u[ind[block_mask]]

    if nrank==2 and rank==1:
        trace_data = None
        block, topRight, bottomLeft, rhsBlock, uBlock = comm.recv(source=0)

    trace, rhsTrace, uTrace = comm.bcast(trace_data, root=0)

    return block, trace, topRight, bottomLeft, rhsBlock, rhsTrace, uBlock, uTrace