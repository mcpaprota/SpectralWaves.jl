# Package Guide

[`SpectralWaves.jl`](https://github.com/mcpaprota/SpectralWaves.jl) is a Julia package for simulation of nonlinear waves propagating over arbitrary topography.

## Instalation

```julia
pkg> add SpectralWaves
julia> using SpectralWaves
```

## Quick start

We begin our introduction with an evolution of a free surface for an initial bump of water in a domain of length ``\ell`` and mean depth ``d`` corresponding to the still water level.

```@example 0
using SpectralWaves
using CairoMakie # plotting package

d = 1.0 # water depth (m)
ℓ = 10 # fluid domain length (m)
nothing # hide
```

We define a number of numerical model parameters.

```@example 0
ℐ = 100 # number of harmonics
Δt = 0.01 # time step (s)
t₀ = 0.0 # initial time (s)
τ = 10.0 # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range
nothing # hide
```

We initialize wave problem using a struct `p::Problem`.

```@example 1
p = Problem(ℓ, d, ℐ, t)
nothing # hide
```

