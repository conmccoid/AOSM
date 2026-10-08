# trace-AOSM -- basic version of multi-subdomain AOSM
# We assume that input is the blocks of the larger matrix, with the trace already identified.

import numpy as np

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
        self.block = block # diagonal block for this subdomain
        self.trace = trace # trace block of the matrix
        self.topRight = topRight # top right blocks of the matrix
        self.bottomLeft = bottomLeft # bottom left blocks of the matrix
        self.rhsBlock = rhsBlock # right hand sides for each block
        self.rhsTrace = rhsTrace # right hand side for the trace
        mod = np.linalg.solve(self.block, self.rhsBlock)
        self.rhsMod = np.dot(self.bottomLeft, mod)

    def setup(self):
        self.T0 = np.zeros((self.sizeTrace, self.sizeTrace)) # initial T matrix
        self.T = self.T0.copy()
        self.S = np.zeros((self.sizeTrace, self.sizeTrace))
        self.wBlock = []
        self.wTrace = []
        self.V = []

    def formT(self):
        T = self.T0.copy()
        for i in range(len(self.wTrace)):
            T -= np.outer(self.V[i], self.wTrace[i])
        return T
    # this should eventually be replaced with the Schur shuffle or similar

    def updateT(self):
        self.T = self.formT()
        # nb: produces subpar results
    
    def solveBlock(self, T, rhs):
        # lin. solve of a specific block
        # can also be made nonlin?
        matrix = np.block([[self.block, self.topRight],[self.bottomLeft, self.trace + T]])
        return np.linalg.solve(matrix, rhs)

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
        rhsModTotal = np.zeros(self.sizeTrace)
        comm.Allreduce(self.rhsMod, rhsModTotal, op=MPI.SUM)
        ftri = self.rhsTrace - rhsModTotal + 2*self.rhsMod

        Wit = self.wTrace
        Wii = self.wBlock
        Vi = self.V
        n = self.sizeBlock

        # initial guesses
        uit = ftri
        uii = -np.linalg.solve(Aii, np.dot(Ait, uit))
        f = np.concatenate((fi, ftri - np.dot(Ati, uii) + np.dot(self.T0, uit)))
        ui = self.solveBlock(self.formT(), f)

        # initial difference vectors
        dii = ui[:n] - uii
        dit = ui[n:] - uit
        a = np.linalg.norm(dit)
        Wit.append(dit/a)
        Wii.append(dii/a)
        Vi.append(-np.dot(Ati, Wii[0]) + np.dot(self.T0, Wit[0]))

        # main loop for MGS orthogonalization
        while a > tol and len(Wit) < maxit:
            di = self.solveBlock(self.formT(), np.concatenate((np.zeros(n), a*Vi[-1]))) # nb: need to replace formT for accurate and fast computation
            ui += di # update solution
            dii = di[:n]
            dit = di[n:] # difference on trace
            a = self.MGS(Wit, Wii, Vi, dit, dii, Ati) # MGS orthogonalization

        uBlock = ui[:n]
        uTrace = ui[n:]
        return uBlock, uTrace, len(Wit), a # return updated solution for the subdomain

    def constructS(self, blockIndex):
        for i in range(self.nBlocks):
            if i != blockIndex:
                self.S[blockIndex] += self.T[i]
    # nb: no longer necessary?
    # MPI: if it is necessary, needs complete overhaul

    def formResidual(self, uBlock, uTrace):
        rtrTotal = -np.matmul(self.bottomLeft, uBlock)
        rtrTrace = -np.matmul(self.trace, uTrace)
        comm.Allreduce(MPI.IN_PLACE, rtrTotal, op=MPI.SUM)
        comm.Allreduce(MPI.IN_PLACE, rtrTrace, op=MPI.SUM)
        rtr = self.rhsTrace + rtrTotal + rtrTrace/size
        return rtr

    def main(self,tol,maxit):
        # initial guess
        # problem parameters

        self.setup() # setup matrices and vectors

        # prepare for global iteration
        uBlock, uTrace, N_it, _ = self.adaptTransmission(tol, maxit) # find adapted transmission conditions
        self.updateT() # update the T matrix for this block
        N_total = comm.allreduce(N_it, op=MPI.SUM) # sum number of iterations from all processes
        if rank==0:
            print(f"Total number of iterations to adapt TCs: {N_total}")
        comm.Barrier() # synchronize processes

        # gather all T matrices from all processes
        T = np.zeros((self.sizeTrace, self.sizeTrace)) # initialize T matrix for summation
        comm.Allreduce(self.T, T, op=MPI.SUM) # gather T matrices from all processes
        self.S = T - self.T # update S matrices for each block

        rtr = self.formResidual(uBlock, uTrace)

        # main loop
        counter=0
        while np.linalg.norm(rtr) > 1e-8: # while residual is large
            # build right hand side
            rhsTraceTotal = np.zeros(self.sizeTrace)
            rhsTraceTotal = np.matmul(self.T, uTrace) - np.matmul(self.bottomLeft, uBlock)
            comm.Allreduce(MPI.IN_PLACE, rhsTraceTotal, op=MPI.SUM)
            rhsTraceTotal += self.rhsTrace
            rhsTrace = rhsTraceTotal - (np.matmul(self.T, uTrace) - np.matmul(self.bottomLeft, uBlock))

            u = self.solveBlock(self.S, np.concatenate((self.rhsBlock, rhsTrace))) # solve the global system
            uBlock = u[:self.sizeBlock] # store new solutions
            uTrace = u[self.sizeBlock:]

            rtr = self.formResidual(uBlock, uTrace) # update residual
            counter+=1 # increment counter
            if rank==0:
                print(f"Iteration: {counter}, Residual: {np.linalg.norm(rtr)}")

        # generate the final trace solution by averaging across processes
        comm.Allreduce(MPI.IN_PLACE, uTrace, op=MPI.SUM) # gather the averaged trace solution across processes
        uTrace /= size # average the trace solution across all processes
        return uBlock, uTrace