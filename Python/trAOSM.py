# trace-AOSM -- basic version of multi-subdomain AOSM
# We assume that input is the blocks of the larger matrix, with the trace already identified.

import numpy as np

# for an eventual sparse version
from scipy.sparse.linalg import spsolve
from scipy import sparse

class AOSM:
    def __init__(self, blocks, trace, topRight, bottomLeft, rhsBlocks, rhsTrace):
        self.nBlocks = len(blocks) # number of blocks
        self.sizeTrace = np.shape(trace)[0] # size of trace (might be wrong size)
        self.blocks = blocks # diagonal blocks of the matrix
        self.trace = trace # trace block of the matrix
        self.topRight = topRight # top right blocks of the matrix
        self.bottomLeft = bottomLeft # bottom left blocks of the matrix
        self.rhsBlocks = rhsBlocks # right hand sides for each block
        self.rhsTrace = rhsTrace # right hand side for the trace
        self.rhsMods = []
        self.sizeBlocks = []
        for i in range(self.nBlocks):
            self.sizeBlocks.append(np.shape(self.blocks[i])[0]) # size of each block
            mod1 = np.linalg.solve(self.blocks[i], self.rhsBlocks[i])
            mod2 = np.dot(self.bottomLeft[i], mod1)
            self.rhsMods.append(mod2) # modifications to rhsTrace

    def setup(self):
        self.T0 = np.zeros((self.sizeTrace, self.sizeTrace)) # initial T matrix
        self.T = [self.T0.copy() for _ in range(self.nBlocks)]
        self.S = [np.zeros((self.sizeTrace, self.sizeTrace)) for _ in range(self.nBlocks)]
        self.wBlocks = [[] for _ in range(self.nBlocks)]
        self.wTrace = [[] for _ in range(self.nBlocks)]
        self.V = [[] for _ in range(self.nBlocks)]

    def formT(self, blockIndex):
        self.T[blockIndex] = self.T0.copy()
        for i in range(len(self.wTrace[blockIndex])):
            self.T[blockIndex] -= np.outer(self.V[blockIndex][i], self.wTrace[blockIndex][i])
    # this should eventually be replaced with the Schur shuffle or similar

    def updateT(self, blockIndex, w, v):
        self.T[blockIndex] -= np.outer(v, w)
        # nb: produces subpar results
    
    def solveBlock(self, blockIndex, T, rhs):
        # lin. solve of a specific block
        # can also be made nonlin?
        matrix = np.block([[self.blocks[blockIndex], self.topRight[blockIndex]], [self.bottomLeft[blockIndex], self.trace + T]])
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

    def adaptTransmission(self, blockIndex, tol, maxit):
        # adapt transmission conditions for a specific block

        # rename variables for brevity
        Aii = self.blocks[blockIndex]
        Ait = self.topRight[blockIndex]
        Ati = self.bottomLeft[blockIndex]
        Ti  = self.T[blockIndex]
        fi = self.rhsBlocks[blockIndex]
        ftri = self.rhsTrace + self.rhsMods[blockIndex] - np.sum([self.rhsMods[i] for i in range(self.nBlocks) if i != blockIndex], axis=0) # replaces modRhs, which was only called at this line
        Wit = self.wTrace[blockIndex]
        Wii = self.wBlocks[blockIndex]
        Vi = self.V[blockIndex]
        n = self.sizeBlocks[blockIndex]

        # initial guesses
        uit = ftri
        uii = -np.linalg.solve(Aii, np.dot(Ait, uit))
        f = np.concatenate((fi, ftri - np.dot(Ati, uii) + np.dot(Ti, uit)))
        ui = self.solveBlock(blockIndex, Ti, f)

        # initial difference vectors
        dii = ui[:n] - uii
        dit = ui[n:] - uit
        a = np.linalg.norm(dit)
        Wit.append(dit/a)
        Wii.append(dii/a)
        Vi.append(-np.dot(Ati, Wii[0]) + np.dot(Ti, Wit[0]))

        # main loop for MGS orthogonalization
        while a > tol and len(Wit) < maxit:
            di = self.solveBlock(blockIndex, Ti, np.concatenate((np.zeros(n), a*Vi[-1]))) # nb: need to replace formT for accurate and fast computation
            ui += di # update solution
            dii = di[:n]
            dit = di[n:] # difference on trace
            a = self.MGS(Wit, Wii, Vi, dit, dii, Ati) # MGS orthogonalization
            self.formT(blockIndex) # update T matrix for this block

        uBlock = ui[:n]
        uTrace = ui[n:]
        return uBlock, uTrace, len(Wit), a # return updated solution for the subdomain

    def constructS(self, blockIndex):
        for i in range(self.nBlocks):
            if i != blockIndex:
                self.S[blockIndex] += self.T[i]

    def formResidual(self, uBlocks, uTrace):
        rtr = self.rhsTrace.copy()
        for i in range(self.nBlocks):
            rtr -= np.matmul(self.bottomLeft[i], uBlocks[i])
        rtr -= np.matmul(self.trace, uTrace)
        return rtr

    # ready for global iteration?
    def solveGlobal(self, uBlocks, uTraces):
        uTrace = np.mean(uTraces, axis=0) # average the trace solutions
        rtr = self.formResidual(uBlocks, uTrace)

        # main loop
        counter=0
        while np.linalg.norm(rtr) > 1e-8: # while residual is large
            # initialize storage for new solutions
            uBlocks_new = [[] for _ in range(self.nBlocks)]
            uTrace_new = [[] for _ in range(self.nBlocks)]
            for i in range(self.nBlocks): # for every subdomain
                # build the right hand side
                rhsTrace = self.rhsTrace.copy()
                for j in range(self.nBlocks):
                    if j != i:
                        Ti = self.T[j]
                        rhsTrace += np.matmul(Ti, uTraces[j]) - np.matmul(self.bottomLeft[j], uBlocks[j])
                ui = self.solveBlock(i, self.S[i], np.concatenate((self.rhsBlocks[i],rhsTrace))) # solve the subdomain
                uBlocks_new[i] = ui[:self.sizeBlocks[i]] # store new solutions
                uTrace_new[i] = ui[self.sizeBlocks[i]:]
            # update solution and residual
            uBlocks = uBlocks_new
            uTraces = uTrace_new
            uTrace = np.mean(uTrace_new, axis=0) # average the trace solutions
            rtr = self.formResidual(uBlocks, uTrace)
            counter+=1
            print(f"Iteration: {counter}, Residual: {np.linalg.norm(rtr)}")
        return uBlocks, uTrace

    def main(self,tol,maxit):
        # initial guess
        # problem parameters

        self.setup() # setup matrices and vectors

        uBlocks = []
        uTrace = []
        N_it = 0.0
        # prepare for global iteration
        for i in range(self.nBlocks):
            uBlock, uTrace_i, n_it, _ = self.adaptTransmission(i, tol, maxit) # find adapted transmission conditions
            uBlocks.append(uBlock) # store the temporary solution for each block
            uTrace.append(uTrace_i)
            N_it += n_it
        print(f"Total number of iterations to adapt: {N_it}")

        for i in range(self.nBlocks):
            self.constructS(i) # construct S matrices for each block

        uBlocks, uTrace = self.solveGlobal(uBlocks, uTrace) # solve the global problem

        return uBlocks, uTrace # return the solution