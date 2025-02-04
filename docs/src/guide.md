# Package Guide

[`SpectralWaves.jl`](https://github.com/mcpaprota/SpectralWaves.jl) is a Julia package for simulation of nonlinear waves propagating over arbitrary topography.

## Instalation

```julia
pkg> add SpectralWaves
julia> using SpectralWaves
```

## Quick start

We begin our introduction with an evolution of a free surface ``\eta(x, t)`` for an initial bump of water in a domain of length ``\ell`` and still water depth ``d``.

```@example 0
using SpectralWaves
using CairoMakie # plotting package

d = 1.0 # water depth (m)
ℓ = 10 # fluid domain length (m)
nothing # hide
```

We define a number of numerical model parameters.

```@example 0
ℐ = 20 # number of harmonics
Δt = 0.01 # time step (s)
t₀ = 0.0 # initial time (s)
τ = 2.0 # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range
nothing # hide
```

We initialize wave problem `p` using a struct `Problem`.

```@example 0
p = Problem(ℓ, d, ℐ, t)
nothing # hide
```

The free surface corresponds to a Gaussian `surface_bump!` of characteristic height ``h`` and length ``\lambda``.

```@example 0
h = 0.3d # bump height (m)
λ = 0.2ℓ # bump length (m)
surface_bump!(p, h, λ)
nothing # hide
```

We solve the problem.

```@example 0
solve_problem!(p)
nothing # hide
```

Finally, we may calculate free surface elevation using `water_surface` function for a range of spatial points `x`

```@example 0
x = range(- ℓ / 2, ℓ / 2, length = 1001) # spatial range
η(x, n) = water_surface(p, x, n)
nothing # hide
```
and plot the results.

```@example 0
set_theme!(theme_latexfonts()) # set latex fonts
fig = Figure(size = (600, 300)) # initialize a figure
ax = Axis(fig[1, 1], 
        xlabel = L"$x$ (m)", 
        ylabel = L"$z$ (m)") # define axis with labels
band!(ax, x, η.(x, lastindex(t)), -d, 
        color=:azure) # plot final water bulk
lines!(ax, x, η.(x, lastindex(t)), 
        color=:black, 
        linewidth = 1, 
        label=L"\eta(\tau=2\,\mathrm{s})") # plot final free surface line
band!(ax, x, -1.1d, - d, 
        color=:wheat) # plot bottom bulk
hlines!(- d, 
        color=:black, 
        linewidth = 0.7) # plot bottom line
limits!(ax, x[1], x[end], -1.1d, d) # set limits
fig # display figure
```

Now, we are going to add some bottom variation by applying `bottom_bump!` to a fresh problem p2. Please note that we set bottom nonlinearity parameter `M_b=30`, while we initialize the free surface with our previous settings using `surface_bump!`.

```@example 0
p2 = Problem(ℓ, d, ℐ, t, M_b=30)
surface_bump!(p2, h, λ)
h₀, λ₀, x₀ = 0.9d, 0.5ℓ, 0.5ℓ # bottom bump height, length, and offset
bottom_bump!(p2, h₀, λ₀, x₀)
```

Again, we solve the problem.

```@example 0
solve_problem!(p2)
nothing # hide
```

We calculate free surface elevation and bottom topography using `water_surface` and `bottom_surface` for a range of points

```@example 0
x = range(- ℓ / 2, ℓ / 2, length = 1001) # spatial range
η(x, n) = water_surface(p2, x, n)
β(x) = bottom_surface(p2, x)
nothing # hide
```
and plot the results.

```@example 0
set_theme!(theme_latexfonts()) # set latex fonts
fig = Figure(size = (600, 300)) # initialize a figure
ax = Axis(fig[1, 1], 
        xlabel = L"$x$ (m)", 
        ylabel = L"$z$ (m)") # define axis with labels
band!(ax, x, η.(x, lastindex(t)), -d, 
        color=:azure) # plot final water bulk
lines!(ax, x, η.(x, lastindex(t)), 
        color=:black, 
        linewidth = 1, 
        label=L"\eta(\tau=2\,\mathrm{s})") # plot final free surface line
band!(ax, x, β.(x) .- d, - 1.1d, 
        color=:wheat) # plot bottom bulk
lines!(ax, x, β.(x) .- d, 
        color=:black, 
        linewidth = 1) # plot bottom line
limits!(ax, x[1], x[end], -1.1d, d) # set limits
fig # display figure
```