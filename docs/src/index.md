```@meta
CurrentModule = SpectralWaves
```

# SpectralWaves.jl

A Fourier Galerkin method solution to nonlinear waves propagating over topography.

## Overview

[`SpectralWaves.jl`](https://github.com/mcpaprota/SpectralWaves.jl) is a Julia package for simulation of nonlinear waves propagating over arbitrary topography under potential flow assumptions. The solution is derived using a Fourier Galerkin spectral method in terms of amplitudes of free-surface elevation and velocity potential, while inverse Fourier transform is used to get a phase-resolved wave field. Four wave generation mechanisms are supported - initial conditions, linear wavemaker forcing, pressure forcing (to be implemented), moving bottom.

## Wave problem

We consider waves propagating over arbitrary bottom topography in a periodic fluid domain of length $\ell$ and mean depth $d$ (corresponding to still water level). A Cartesian coordinate system is used to define fluid elements along horizontal $x$-axis coinciding with undisturbed free surface and upward-pointing and vertical $z$-axis. Undulating free surface is described by means of $\eta(x, t)$, while bottom topography is considered as fluctuations $\beta(x, t)$ around $-d$. 


Documentation for [SpectralWaves].
