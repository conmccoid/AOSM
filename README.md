# Adaptive Optimized Schwarz Methods (AOSMs)

This repository collects the latest versions of AOSMs.

Project goals:
- [ ] fully functioning multi-subdomain AOSM
- [ ] crosspoints
- [ ] symmetrized cells
- [ ] alternative T decomps
- [ ] combine with multigrid for scaling
- [ ] high dimension nesting
- [ ] analysis
    - [ ] Krylov type
    - [ ] convergence
    - [ ] efficiency / flop count
    - [ ] Schur complements
    - [ ] dependence on initial guess / choice of initial guess
- [ ] applications
    - [ ] nonlinear (Newton & quasi-Newton)
    - [ ] time dependent (Runge-Kutta & exact time differencing)
    - [ ] machine learning (optimization of T / AOSM speed-ups)

MATLAB prototypes available:
- [x] trAOSM, multi-subdomain full trace decomposition
- [ ] restricted trace
- [ ] solving exclusively on the trace then propagating
- [ ] diagonally dominant decomposition
- [ ] hierarchical decomposition

## Organization

The repository is organized by code base.
- MATLAB: prototyping new methods.
- Python: allows testing of parallelization in more contexts, as well as compatibility with FEniCSx and other FEM software.
- C: allows compatibility with PETSc.

- [x] Decide on additional code bases, such as FEniCS and PETSc

## Notes

Examples with crosspoints fail because the trace has an interior, i.e. the entire stencil of one or more points lies in the trace.
By giving crosspoints their own subdomain, the method converges.
However, it appears that the method is then unstable and heavily reliant on the initial guess surrounding the crosspoints.

Diagonal dominance in ATCs must be taken advantage of.

### Parallelization of code

I need a consistent workflow for parallel computing.
I can't seem to find any good advice on this, it seems like there's no one way to do things.

One option is to create one subdomain per rank, with the trace copied to all ranks.
The algorithm would then operate first on each subdomain to obtain the transmission conditions,
then transmit information between all ranks,
and finally compute a complete solution by iterating globally.

This is distinct from all serial incarnations of the code, which has the subdomains on all ranks and then splits up the computation from there.

The thing I'm mentally stuck on, I think, is what to do when the number of ranks is not equal to the number of subdomains.
In particular, if I want to test the code in serial, the first option above won't be able to run over all of the subdomains (?).

Perhaps the best workflow is to assume the first option, which should be the fastest, but then patch it for serial runs and rank less than number of blocks.