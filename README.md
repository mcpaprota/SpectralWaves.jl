# SpectralWaves

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mcpaprota.github.io/SpectralWaves.jl/dev/)
[![Build Status](https://github.com/mcpaprota/SpectralWaves.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mcpaprota/SpectralWaves.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![DOI](https://zenodo.org/badge/909337001.svg)](https://doi.org/10.5281/zenodo.14998246)

A Julia implementation of a spectral method for nonlinear water waves propagating over topography.

## Install

```julia
pkg> add SpectralWaves
julia> using SpectralWaves
```

## Documentation

The very first version of documentation is available [here](https://mcpaprota.github.io/SpectralWaves.jl/dev/).

## Citation
Please cite as:
Maciej Paprota,
A fully spectral framework for nonlinear water waves propagating over topography,
Coastal Engineering,
Volume 200,
2025,
104759,
ISSN 0378-3839,
https://doi.org/10.1016/j.coastaleng.2025.104759

The paper is freely available until June, 19, 2025 under this [sharelink](https://authors.elsevier.com/a/1l0cv_8cdlyPzW).


## Related packages

Below is a list of packages, which inspired development of `SpectralWaves.jl`:

- [`DSP.jl`](https://github.com/JuliaDSP/DSP.jl) - digital signal processing tools;
- [`FourierFlows.jl`](https://github.com/FourierFlows/FourierFlows.jl) - Fourier-collocation methods for partial differential equations;
- [`WaterWaves1D.jl`](https://github.com/WaterWavesModels/WaterWaves1D.jl) - a collection of one-dimensional wave models;
- [`Oceananigans.jl`](https://github.com/CliMA/Oceananigans.jl) - solvers for ocean dynamics equations.
