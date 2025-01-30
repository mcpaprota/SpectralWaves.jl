# Package Guide

[`SpectralWaves.jl`](https://github.com/mcpaprota/SpectralWaves.jl) is a Julia package for simulation of nonlinear waves propagating over arbitrary topography.

## Instalation

```julia
pkg> add SpectralWaves
julia> using SpectralWaves
```

## Quick start

We begin our introduction with application of the model to linear and regular waves of length ``L`` and height ``H`` propagating in water of constant depth ``d``. We consider one wave along the length of the domain ``\ell``. 

```@example 1
using SpectralWaves
using CairoMakie # plotting package

L = 2.0 # wavelength (m)
H = 0.1 # wave height (m)
d = 1.0 # water depth (m)
ℓ = L # fluid domain length (m) - one wave
nothing # hide
```

We need a wave period ``T``.

```@example 1
k = 2π / L # wave number (rad/m)
ω = sqrt(g * k * tanh(k * d)) # angular wave frequency (rad/s)
T = 2π / ω # wave period (s)
```

We define a number of numerical model parameters.

```@example 1
ℐ = 1 # number of harmonics
nΔt = 200 # number of time increments per wave period
Δt = T / nΔt # time step (s)
nT = 1 # number of periods
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range
nothing # hide
```

We initialize wave problem using a struct `p::Problem`.

```@example 1
p = Problem(ℓ, d, ℐ, t)
nothing # hide
```

Initial condition values of ``\hat{\eta}``, ``\dot{\eta}``, ``\hat{\phi}``, and ``\dot{\phi}`` are computed and inserted into vectors `η̂`, `η̇`, `ϕ̂`, `ϕ̇` using `linear_regular_wave!` in-place function.

```@example 1
linear_regular_wave!(p, H, ω)
nothing # hide
```

Now, we are ready to solve a problem. We use an in-place function `solve_problem!` which stores the values of solution coefficients in vectors `η̂`, `η̇`, `ϕ̂`, `ϕ̇`, `ψ̂`, `ψ̇`. In our case only `η̂` will be further processed.

```@example 1
solve_problem!(p)
nothing # hide
```

Finally we are ready to plot evolution of absolute, real and imaginary values of a complex wave amplitude `η̂` over time `t` corresponding to one wave period.

```@example 1
set_theme!(theme_latexfonts())
fig = Figure(size = (400, 300))
ax = Axis(fig[1, 1], xlabel = L"t/T", ylabel = L"4η̂/H", xticks = 0:0.1:p.N, yticks = -1:0.5:1)
lines!(ax, t / T, 4 * abs.(p.η̂[1, p.O:end]) / H, label = L"|η̂|")
lines!(ax, t / T, 4 * real.(p.η̂[1, p.O:end]) / H, label = L"Re(η̂)")
lines!(ax, t / T, 4 * imag.(p.η̂[1, p.O:end]) / H, label = L"Im(η̂)")
axislegend(ax, position = :lb)
fig
```

Now, if we want to plot a time-series of free-surface elevation at some location (e.g. corresponding to the middle of the fluid domain ``x_0 = \ell/2``) ``\eta(t[n], x_0)``, we use an `inverse_fourier_transform` function, like so

```@example 1
x_0 = ℓ/2
η(n) = inverse_fourier_transform(p.η̂[:, n], p.κ, x_0)
fig = Figure(size = (400, 300))
ax = Axis(fig[1, 1], xlabel = L"$t$ (s)", ylabel = L"$\eta$ (m)")
lines!(ax, t, η.(p.O:p.N+p.O), color = :blue, linewidth = 2)
limits!(ax, 0, T, -H, H)
fig
```
