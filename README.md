# Adaptive Optimized Schwarz Methods (AOSMs)

This repository collects the latest versions of AOSMs.

Project goals:
- [ ] fully functioning multi-subdomain AOSM
- [ ] crosspoints
- [ ] symmetrized cells
- [ ] combine with multigrid for scaling
- [ ] high dimension nesting

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

Diagonal dominance in ATCs must be taken advantage of.