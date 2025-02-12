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

- [ ] Decide on additional code bases, such as FEniCS and PETSc

## Notes

Examples with crosspoints fail because the trace has an interior, i.e. the entire stencil of one or more points lies in the trace.
By giving crosspoints their own subdomain, the method converges.
However, it appears that the method is then unstable and heavily reliant on the initial guess surrounding the crosspoints.

Diagonal dominance in ATCs must be taken advantage of.