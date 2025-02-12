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

Currently can handle 3 subdomains very well, loves strips.
4 subdomains breaks, likely due to crosspoint.
- Idea: crosspoint is a problem because the entire stencil of this point lies in the trace. Why should it be that the trace can't resolve interior points? This is worth investigating.

Diagonal dominance in ATCs must be taken advantage of.