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
    
    def solveBlock(self, blockIndex, rhs):
        # lin. solve of a specific block
        # can also be made nonlin?
        block = self.blocks[blockIndex]
        topRight = self.topRight[blockIndex]
        bottomLeft = self.bottomLeft[blockIndex]
        T = self.formT(blockIndex)
        matrix = np.block([[block, topRight], [bottomLeft, self.trace + T]])
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
        Aii = self.blocks[blockIndex]
        Ati = self.topRight[blockIndex]
        Ait = self.bottomLeft[blockIndex]
        fi = self.rhsBlocks[blockIndex]
        ftri = self.modRhs(blockIndex)
        Wit = self.wTrace[blockIndex]
        Wii = self.wBlocks[blockIndex]
        Vi = self.V[blockIndex]
        n = self.sizeBlocks[blockIndex]

        uit = ftri
        uii = -np.linalg.solve(Aii, np.linalg.matmul(Ait, uit))
        f = np.concatenate((fi, ftri - np.matmul(Ati, uii) + np.matmul(self.T, uit)))
        ui = self.solveBlock(blockIndex, f)

        dii = ui[:n] - uii
        dit = ui[n:] - uit
        a = np.norm(dit)
        Wit.append(dit/a)
        Wii.append(dii/a)
        Vi.append(-np.matmul(Ati, Wii[0]) + np.matmul(self.T, Wit[0]))

        while a > 1e-12 and len(Wit) < 1000:
            di = self.solveBlock(blockIndex, np.concatenate((np.zeros(n), a*Vi[-1])))
            dii = di[:n]
            dit = di[n:]
            a = self.MGS(Wit, Wii, Vi, dit, dii, Ati)
        # this version currently does not adapt the solution, only the transmission conditions

    def constructS(self, blockIndex):
        S = self.S[blockIndex]
        for i in range(self.nBlocks):
            if i != blockIndex:
                S += self.formT(i)

    # ready for global iteration?