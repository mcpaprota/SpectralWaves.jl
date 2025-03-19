# SpectralWaves

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mcpaprota.github.io/SpectralWaves.jl/dev/)
[![Build Status](https://github.com/mcpaprota/SpectralWaves.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mcpaprota/SpectralWaves.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![DOI](https://zenodo.org/badge/909337001.svg)](https://doi.org/10.5281/zenodo.14998246)

A Julia implementation of a spectral method for nonlinear water waves propagating over topography.

> [!WARNING]  
> There is a bug somewhere in the code. It makes some results of nonlinear wave problems inaccurate. We are working on it.

## Install

```julia
pkg> add SpectralWaves
julia> using SpectralWaves
```

## Documentation

The very first version of documentation is available [here](https://mcpaprota.github.io/SpectralWaves.jl/dev/).

## Related packages

Below is a list of packages, which inspired development of `SpectralWaves.jl`:

- [`DSP.jl`](https://github.com/JuliaDSP/DSP.jl) - digital signal processing tools;
- [`FourierFlows.jl`](https://github.com/FourierFlows/FourierFlows.jl) - Fourier-collocation methods for partial differential equations;
- [`WaterWaves1D.jl`](https://github.com/WaterWavesModels/WaterWaves1D.jl) - a collection of one-dimensional wave models;
- [`Oceananigans.jl`](https://github.com/CliMA/Oceananigans.jl) - solvers for ocean dynamics equations.
