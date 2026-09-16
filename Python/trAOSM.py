# trace-AOSM -- basic version of multi-subdomain AOSM
# We assume that input is the blocks of the larger matrix, with the trace already identified.

import numpy as np

class AOSM:
    def __init__(self, blocks, trace, topRight, bottomLeft, rhsBlocks, rhsTrace):
        self.nBlocks = len(blocks) # number of blocks
        self.sizeTrace = np.size(trace)[0] # size of trace (might be wrong size)
        self.blocks = blocks # diagonal blocks of the matrix
        self.trace = trace # trace block of the matrix
        self.topRight = topRight # top right blocks of the matrix
        self.bottomLeft = bottomLeft # bottom left blocks of the matrix
        self.rhsBlocks = rhsBlocks # right hand sides for each block
        self.rhsTrace = rhsTrace # right hand side for the trace
        self.rhsMods = []
        self.sizeBlocks = []
        for i in range(self.nBlocks):
            self.sizeBlocks.append(np.size(self.blocks[i])[0]) # size of each block
            mod1 = np.linalg.solve(self.blocks[i], self.rhsBlocks[i])
            mod2 = np.linalg.matmul(self.bottomLeft[i], mod1)
            self.rhsMods.append(mod2) # modifications to rhsTrace

    def setup(self):
        self.T = np.empty(self.sizeTrace,self.sizeTrace)
        self.wBlocks = [[] for _ in range(self.nBlocks)]
        self.wTrace = [[] for _ in range(self.nBlocks)]
        self.V = [[] for _ in range(self.nBlocks)]
        self.S = [[] for _ in range(self.nBlocks)]

    def modRhs(self, blockIndex):
        return self.rhsTrace + self.rhsMods[blockIndex] - np.sum([self.rhsMods[i] for i in range(self.nBlocks) if i != blockIndex], axis=0)

    def formT(self, blockIndex):
        T = self.T
        for i in range(len(self.wTrace[blockIndex])):
            T -= np.outer(self.V[blockIndex][i], self.wTrace[blockIndex][i])
        return T
    # this should eventually be replaced with the Schur shuffle or similar
    
    def solveBlock(self, blockIndex, T, rhs):
        # lin. solve of a specific block
        # can also be made nonlin?
        matrix = np.block([[self.blocks[blockIndex], self.topRight[blockIndex]], [self.bottomLeft[blockIndex], self.trace + T]])
        return np.linalg.solve(matrix, rhs)

    def MGS(self, Wit, Wii, Vi, dit, dii, Ati):
        # Modified Gram-Schmidt orthogonalization
        Wii.append(dii)
        Wit.append(dit)
        Vi.append(-np.linalg.matmul(Ati, dii) + np.linalg.matmul(self.T,dit))
        for k in range(len(Wit)-1):
            r = np.dot(Wit[k], Wit[-1])
            Wit[-1] -= r * Wit[k]
            Wii[-1] -= r * Wii[k]
            Vi[-1] -= r * Vi[k]
        a = np.linalg.norm(Wit[-1])
        if a > 1e-12:
            Wit[-1] /= a
            Wii[-1] /= a
            Vi[-1] /= a
        return a
        # not sure if this function will correctly modify Wit, Wii and Vi externally

    def adaptTransmission(self, blockIndex):
        # adapt transmission conditions for a specific block

        # rename variables for brevity
        Aii = self.blocks[blockIndex]
        Ati = self.topRight[blockIndex]
        Ait = self.bottomLeft[blockIndex]
        fi = self.rhsBlocks[blockIndex]
        ftri = self.modRhs(blockIndex)
        Wit = self.wTrace[blockIndex]
        Wii = self.wBlocks[blockIndex]
        Vi = self.V[blockIndex]
        n = self.sizeBlocks[blockIndex]

        # initial guesses
        uit = ftri
        uii = -np.linalg.solve(Aii, np.linalg.matmul(Ait, uit))
        f = np.concatenate((fi, ftri - np.matmul(Ati, uii) + np.matmul(self.T, uit)))
        ui = self.solveBlock(blockIndex, self.formT(blockIndex), f)

        # initial difference vectors
        dii = ui[:n] - uii
        dit = ui[n:] - uit
        a = np.norm(dit)
        Wit.append(dit/a)
        Wii.append(dii/a)
        Vi.append(-np.matmul(Ati, Wii[0]) + np.matmul(self.T, Wit[0]))

        # main loop for MGS orthogonalization
        while a > 1e-12 and len(Wit) < 1000:
            di = self.solveBlock(blockIndex, self.formT(blockIndex),np.concatenate((np.zeros(n), a*Vi[-1])))
            dii = di[:n]
            dit = di[n:]
            a = self.MGS(Wit, Wii, Vi, dit, dii, Ati)
        # this version currently does not adapt the solution, only the transmission conditions

    def constructS(self, blockIndex):
        for i in range(self.nBlocks):
            if i != blockIndex:
                self.S[blockIndex] += self.formT(i)

    def formResidual(self, uBlocks, uTrace):
        rtr = self.rhsTrace
        for i in range(self.nBlocks):
            rtr -= np.matmul(self.bottomLeft[i], uBlocks[i])
        rtr -= np.matmul(self.trace, uTrace)
        return rtr

    # ready for global iteration?
    def solveGlobal(self, uBlocks, uTrace):
        rtr = self.formResidual(uBlocks, uTrace)

        # initialize storage for new solutions
        uBlocks_new = uBlocks.copy()
        uTrace_new = uTrace.copy()

        # main loop
        while np.linalg.norm(rtr) > 1e-12: # while residual is large
            for i in range(self.nBlocks): # for every subdomain
                # build the right hand side
                rhs = self.rhsTrace
                for j in range(self.nBlocks):
                    if j != i:
                        Ti = self.formT(j)
                        rhs += np.matmul(Ti, uTrace) - np.matmul(self.bottomLeft[j], uBlocks[j])
                ui = self.solveBlock(i, self.S[i], rhs) # solve the subdomain
                uBlocks_new[i] = ui[:self.sizeBlocks[i]] # store new solutions
                uTrace_new[i] = ui[self.sizeBlocks[i]:]
            # update solution and residual
            uBlocks = uBlocks_new
            uTrace = np.mean(uTrace_new, axis=0) # average the trace solutions
            rtr = self.formResidual(uBlocks, uTrace)
        return uBlocks, uTrace

    def main(self):
        # initial guess
        # problem parameters

        self.setup() # setup matrices and vectors
        
        # prepare for global iteration
        for i in range(self.nBlocks):
            self.adaptTransmission(i) # find adapted transmission conditions
        for i in range(self.nBlocks):
            self.constructS(i) # construct S matrices for each block

        uBlocks, uTrace = self.solveGlobal(uBlocks, uTrace) # solve the global problem