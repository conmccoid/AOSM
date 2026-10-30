# trace-AOSM -- basic version of multi-subdomain AOSM
# We assume that input is the blocks of the larger matrix, with the trace already identified.

import numpy as np
from scipy.linalg import lu_factor, lu_solve

# for an eventual sparse version
from scipy.sparse.linalg import spsolve
from scipy import sparse

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank, size = comm.Get_rank(), comm.Get_size()

class AOSM:
    def __init__(self, block, trace, topRight, bottomLeft, rhsBlock, rhsTrace):
        self.sizeTrace = np.shape(trace)[0] # size of trace (might be wrong size)
        self.sizeBlock = np.shape(block)[0] # size of block

        self.traceRows = np.flatnonzero(np.any(bottomLeft != 0, axis=1))
        self.traceCols = np.flatnonzero(np.any(topRight != 0, axis=0))

        self.block = block # diagonal block for this subdomain
        self.topRight = topRight[:,self.traceCols] # top right blocks of the matrix
        self.bottomLeft = bottomLeft[self.traceRows, :] # bottom left blocks of the matrix
        self.trace = trace[np.ix_(self.traceRows, self.traceCols)] # restricted trace block of the matrix
        self.traceGlobal = trace # global trace block of the matrix

        self.rhsBlock = rhsBlock # right hand sides for each block
        self.rhsTrace = rhsTrace # right hand side for the trace
        mod = np.linalg.solve(self.block, self.rhsBlock)
        self.rhsMod = np.dot(self.bottomLeft, mod)

        allRows = comm.allgather(self.traceRows) # gather row and col indices
        allCols = comm.allgather(self.traceCols)
        self.nbr = {} # initialize list of neighbours (?)
        for j in range(size):
            if j == rank: continue # skip local objects
            # intersect traceRows and allRows[j]
            r, ri, _ = np.intersect1d(self.traceRows, allRows[j], assume_unique=True, return_indices=True)
            c, ci, _ = np.intersect1d(self.traceCols, allCols[j], assume_unique=True, return_indices=True)
            if r.size and c.size:
                # store global indices that are shared between rank 'rank' and rank 'j'
                self.nbr[j] = (ri, ci)

    def setup(self):
        self.T0 = np.zeros(self.trace.shape) # initial T matrix
        self.T = self.T0.copy()
        self.S = np.zeros(self.trace.shape) # initial S matrix
        self.wBlock = []
        self.wTrace = []
        self.V = []
        self.zBlock = []
        self.zTrace = []

    def formT(self):
        T = self.T0.copy()
        for i in range(len(self.wTrace)):
            T -= np.outer(self.V[i], self.wTrace[i])
        return T
    # this should eventually be replaced with the Schur shuffle or similar

    def updateT(self):
        self.T = self.formT()
        # nb: produces subpar results

    def factorBlock(self, T):
        matrix = np.block([[self.block, self.topRight],[self.bottomLeft, self.trace + T]])
        return lu_factor(matrix)

    def MGS(self, Wit, Wii, Vi, dit, dii, Ati):
        # Modified Gram-Schmidt orthogonalization (one step)
        Wii.append(dii)
        Wit.append(dit)
        Vi.append(-np.dot(Ati, dii) + np.dot(self.T0, dit))
        for k in range(len(Wit)-1):
            r = np.dot(Wit[k], Wit[-1])
            Wit[-1] -= r * Wit[k]
            Wii[-1] -= r * Wii[k]
            Vi[-1] -= r * Vi[k]
        a = np.linalg.norm(Wit[-1])
        Wit[-1] /= a
        Wii[-1] /= a
        Vi[-1] /= a
        return a
        # not sure if this function will correctly modify Wit, Wii and Vi externally

    def adaptTransmission(self, tol, maxit):
        # adapt transmission conditions for a specific block

        # rename variables for brevity
        Aii = self.block
        Ait = self.topRight
        Ati = self.bottomLeft
        fi = self.rhsBlock
        rhsModLocal = np.zeros(self.sizeTrace)
        rhsModLocal[self.traceRows] = self.rhsMod
        rhsModTotal = np.zeros(self.sizeTrace)
        comm.Allreduce(rhsModLocal, rhsModTotal, op=MPI.SUM)
        ftri = (self.rhsTrace - rhsModTotal + 2*rhsModLocal)[self.traceRows] # nb: reconsider this...

        Wit = self.wTrace
        Wii = self.wBlock
        Vi = self.V
        Zii = self.zBlock
        Zit = self.zTrace
        n = self.sizeBlock

        # initial guesses
        # uit = ftri
        uit = np.random.rand(self.traceRows.size) # random initial guess for the trace
        uii = -np.linalg.solve(Aii, np.dot(Ait, uit))
        f = np.concatenate((fi, ftri - np.dot(Ati, uii) + np.dot(self.T0, uit)))
        lu = self.factorBlock(self.T0)
        ui = lu_solve(lu, f) # initial solution for the block

        # initial difference vectors
        dii = ui[:n] - uii
        dit = ui[n:] - uit
        a = np.linalg.norm(dit)
        Wit.append(dit/a)
        Wii.append(dii/a)
        Vi.append(-np.dot(Ati, Wii[0]) + np.dot(self.T0, Wit[0]))

        # main loop for MGS orthogonalization
        while a > tol and len(Wit) < maxit:
            zi = lu_solve(lu, np.concatenate((np.zeros(n), Vi[-1]))) # solve main block
            Zii.append(zi[:n])
            Zit.append(zi[n:])
            k = len(Wit)
            Zt = np.column_stack(Zit)
            Z = np.vstack((np.column_stack(Zii), Zt))
            di = Z @ np.linalg.solve(np.eye(k) - np.column_stack(Wit).T @ Zt, np.eye(k)[:,-1]) # Woodbury update
            ui += a*di # update solution
            dii = di[:n]
            dit = di[n:] # difference on trace
            a = a * self.MGS(Wit, Wii, Vi, dit, dii, Ati) # MGS orthogonalization
            print(f"Process {rank}: Iteration {len(Wit)}, norm of difference: {a}", flush=True)

        uBlock = ui[:n]
        uTrace = ui[n:]
        return uBlock, uTrace, len(Wit), a # return updated solution for the subdomain

    def formResidual(self, uBlock, uTrace):
        # nb: self.trace portion should also include all columns which have non-zero value for traceRows
        # rtr = -np.matmul(self.trace,uTrace)
        # send = [None]*size
        # for j, (ri, ci) in self.nbr.items():
        #     send[j] = rtr[ri]
        # recv = comm.alltoall(send)
        # for j, (ri, ci) in self.nbr.items():
        #     rtr[ri] += recv[j]
        # rtr += self.rhsTrace[self.traceRows] - np.matmul(self.bottomLeft, uBlock)

        # residual needs to be global so that subprocesses don't stop until everyone is done
        # otherwise, collectives like comm.alltoall() will hang
        uTraceGlobal = np.zeros((self.sizeTrace,))
        uTraceGlobal[self.traceRows] = uTrace/self.multiplicity
        resTrace = self.traceGlobal @ uTraceGlobal
        resTrace[self.traceRows] += self.bottomLeft @ uBlock
        comm.Allreduce(MPI.IN_PLACE, resTrace, op=MPI.SUM)
        resTrace -= self.rhsTrace

        resBlock = self.block @ uBlock + self.topRight @ uTrace - self.rhsBlock
        res = np.linalg.norm(resBlock)**2
        comm.allreduce(res, op=MPI.SUM) # copilot thinks this discards res?
        res += np.linalg.norm(resTrace)**2
        return np.sqrt(res)

    def sumNeighbours(self, obj, type='row'):
        # sum an object over all neighbouring ranks
        send = [None]*size # initialize data to send to other ranks
        for j, (ri, ci) in self.nbr.items(): # for each neighbouring rank
            if type=='mat':
                send[j] = obj[np.ix_(ri,ci)] # send those parts of obj that are shared
            elif type=='row':
                send[j] = obj[ri]
            elif type=='col':
                send[j] = obj[ci]
        recv = comm.alltoall(send) # send data to appropriate ranks
        out = np.zeros(obj.shape) # initialize object to output results
        for j, (ri, ci) in self.nbr.items():
            if type=='mat':
                out[np.ix_(ri,ci)] += recv[j] # sum up all values from neighbours
            elif type=='row':
                out[ri] += recv[j]
            elif type=='col':
                out[ci] += recv[j]
        return out

    def main(self,tol,maxit):
        # initial guess
        # problem parameters

        self.setup() # setup matrices and vectors

        # prepare for global iteration
        uBlock, uTrace, N_it, _ = self.adaptTransmission(tol, maxit) # find adapted transmission conditions
        self.updateT() # update the T matrix for this block
        N_total = comm.allreduce(N_it, op=MPI.SUM) # sum number of iterations from all processes
        if rank==0:
            print(f"Total number of iterations to adapt TCs: {N_total}", flush=True)
        comm.Barrier() # synchronize processes

        # gather all T matrices from all processes
        S = self.sumNeighbours(self.T, 'mat')
        lu = self.factorBlock(S)

        self.multiplicity = 1 + self.sumNeighbours(np.ones(self.traceCols.size),'col')
        traceCols = self.traceGlobal[:,self.traceCols]

        res = self.formResidual(uBlock, uTrace)

        # main loop
        counter=0
        while res > 1e-8: # while residual is large
            # nb: TRACE IS NO LONGER RESTRICTED IN GLOBAL SOLVE (needs fixing)

            # build right hand side
            # off diagonal portions of Att
            send = [None]*size
            uTraceLocal = uTrace / self.multiplicity # weighted average for prolongation (important when multiplying by Att)
            for j, (ri,ci) in self.nbr.items():
                uTraceGlobal = uTraceLocal.copy()
                uTraceGlobal[ci] = 0 # remove any part associated with i-th trace (I - Ri'Ri)
                send[j] = traceCols @ uTraceGlobal # multiply all rows of Att by neighbouring uTrace
            recv = comm.alltoall(send) # receive neighbouring uTrace contributions
            rhsOffDiagonal = np.zeros(uTrace.shape)
            for j in self.nbr:
                rhsOffDiagonal += recv[j][self.traceRows] # sum up all neighbouring uTrace contributions

            # local transmission conditions
            # temp = self.bottomLeft @ uBlock
            send = [None]*size
            for j, (ri,ci) in self.nbr.items():
                send[j] = self.T[np.ix_(ri,ci)] @ uTrace[ci] - self.bottomLeft[ri,:] @ uBlock
            recv = comm.alltoall(send)
            rhsTrace = np.zeros(uTrace.shape)
            for j, (ri,ci) in self.nbr.items():
                rhsTrace[ri] += recv[j]
            rhsTrace += self.rhsTrace[self.traceRows]
            rhsTrace -= rhsOffDiagonal

            u = lu_solve(lu, np.concatenate((self.rhsBlock, rhsTrace))) # solve the (restricted) global system
            uBlock = u[:self.sizeBlock] # store new solutions
            uTrace = u[self.sizeBlock:] # store new solutions on trace

            res = self.formResidual(uBlock, uTrace) # update residual
            counter+=1 # increment counter
            if rank==0:
                print(f"Iteration: {counter}, Residual: {res}", flush=True)

        # average uTrace over all neighbouring processes
        comm.Barrier()
        uTraceGlobal = np.zeros(self.sizeTrace)
        uTraceGlobal[self.traceCols] = uTrace / self.multiplicity
        comm.Allreduce(MPI.IN_PLACE, uTraceGlobal, op=MPI.SUM)
        return uBlock, uTraceGlobal