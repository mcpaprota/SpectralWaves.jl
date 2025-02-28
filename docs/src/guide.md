# Package Guide

[`SpectralWaves.jl`](https://github.com/mcpaprota/SpectralWaves.jl) is a Julia package for simulation of nonlinear waves propagating over arbitrary topography.

## Instalation

```julia
pkg> add https://github.com/mcpaprota/SpectralWaves.jl
julia> using SpectralWaves
```

## Quick start

We begin our introduction with an evolution of a free surface for an initial bump of water in a domain of length `ℓ` and still water depth `d`. We aim to compare free-surface evolution for the case of constant and uneven bottom.

```@example 0
using SpectralWaves
using CairoMakie # plotting package

d = 1.0 # water depth (m)
ℓ = 10.0 # fluid domain length (m)
nothing # hide
```

We define a number of numerical model parameters.

```@example 0
ℐ = 40 # number of harmonics
Δt = 0.01 # time step (s)
t₀ = 0.0 # initial time (s)
τ = 2.0 # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range
nothing # hide
```

We initialize a constant bottom wave problem `p₀` and an uneven bottom wave problem `p₁` using struct [`Problem`](@ref). Please note that we set a bottom nonlinearity parameter `M_b=40` in case of an uneven bottom, while for constant bottom we leave its default (`M_b=0`) value.

```@example 0
p₀ = Problem(ℓ, d, ℐ, t)
p₁ = Problem(ℓ, d, ℐ, t; M_b=40)
nothing # hide
```

The free surface corresponds to a Gaussian [`surface_bump!`](@ref) of characteristic height `h` and length `λ` and is applied to both problems `p₀` and `p₁`, while we add some bottom variation by applying a Gaussian [`bottom_bump!`](@ref) of characteristic height `h₁` and length `λ₁` to problem `p₁`.

```@example 0
h = 0.4d # surface bump height (m)
λ = 0.1ℓ # surface bump length (m)
surface_bump!(p₀, h, λ)
surface_bump!(p₁, h, λ)
h₁ = 0.9d # bottom bump height (m)
λ₁ = 0.5ℓ # bottom bump length (m)
bottom_bump!(p₁, h₁, λ₁)
nothing # hide
```

We solve both problems.

```@example 0
solve_problem!(p₀)
solve_problem!(p₁)
nothing # hide
```

Finally, we may calculate free surface elevation and bottom surface position using [`water_surface`](@ref) and [`bottom_surface`](@ref) functions

```@example 0
η₀(x, n) = water_surface(p₀, x, n)
η₁(x, n) = water_surface(p₁, x, n)
β(x) = bottom_surface(p₁, x)
nothing # hide
```
and plot the results for a range of spatial points `x`.

```@example 0
x = range(start = - ℓ / 2, stop = ℓ / 2, length = 1001) # spatial range
o₀ = Observable(η₀.(x, firstindex(t))) # set free-surface observable for p₀
o₁ = Observable(η₁.(x, firstindex(t))) # set free-surface observable for p₁
title = Observable(L"t = %$(round(t[1], digits=1))\,\mathrm{s}") # set string observable
set_theme!(theme_latexfonts()) # set latex fonts
fig = Figure(size = (700, 300)) # initialize a figure

# left plot p₀
ax0 = Axis(fig[1, 1], 
        xlabel = L"$x$ (m)", 
        ylabel = L"$z$ (m)") # define axis with labels
band!(ax0, x, o₀, -d, 
        color=:azure) # plot water bulk
lines!(ax0, x, o₀, 
        color=:black, 
        linewidth = 1) # plot free surface line
band!(ax0, x, -1.1d, - d, 
        color=:wheat) # plot bottom bulk
hlines!(ax0, -d, 
        color=:black, 
        linewidth = 0.7) # plot bottom line
limits!(ax0, x[1], x[end], -1.1d, d) # set limits

# right plot p₁
ax1 = Axis(fig[1, 2], 
        xlabel = L"$x$ (m)") # define axis with labels
band!(ax1, x, o₁, β.(x) .-d, 
        color=:azure) # plot water bulk
lines!(ax1, x, o₁, 
        color=:black, 
        linewidth = 1) # plot free surface line
band!(ax1, x, β.(x) .- d, - 1.1d, 
        color=:wheat) # plot bottom bulk
lines!(ax1, x, β.(x) .- d, 
        color=:black, 
        linewidth = 1) # plot bottom line
limits!(ax1, x[1], x[end], -1.1d, d) # set limits
Label(fig[0, :], text=title)

# animate free surface
record(fig, "animation.mp4", 1:lastindex(t);
        framerate = 30) do n
    o₀[] = η₀.(x, n)
    o₁[] = η₁.(x, n)
    title[] = L"t = %$(round(t[n], digits=1))\,\mathrm{s}"
end
nothing # hide
```

```@raw html
<video width="auto" controls autoplay loop>
<source src="../animation.mp4" type="video/mp4">
</video>
```

Please see the next section for more examples.