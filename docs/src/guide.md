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

d = 1.0 # water depth (m)
H = 0.1 # wave height (m)
L = 2.0 # wavelength (m)
ℓ = L # fluid domain length (m) - one wave
nothing # hide
```

We define a number of numerical model parameters.

```@example 1
M_s = 0 # FSBC Taylor series order (linear wave)
M_b = 0 # BBC Taylor series order (horizontal bottom)
ℐ = 1 # number of harmonics
nΔt = 200 # number of time increments per wave period
nT = 1 # number of periods
N = nΔt * nT # number of time steps
nothing # hide
```

We initialize wave problem using `init_problem` function.

```@example 1
κ, η̂, η̇, β̂, β̇, ϕ̂, ϕ̇, ψ̂, ψ̇, p̂, χ, ξ, ζ, 𝒯, 𝒮, O = init_problem(ℓ, d, ℐ, N)
nothing # hide
```

Initial condition values of ``\hat{\eta}``, ``\dot{\eta}``, ``\hat{\phi}``, and ``\dot{\phi}`` are computed and inserted into vectors `η̂`, `η̇`, `ϕ̂`, `ϕ̇` using `linear_regular_wave!` in-place function, which additionally returns values of a wave period ``T`` and a time increment ``\Delta t``.

```@example 1
T, Δt = linear_regular_wave!(η̂, η̇, ϕ̂, ϕ̇, H, L, d, ℐ, nΔt, O)
nothing # hide
```

Now, we are ready to solve a problem. We use an in-place function `solve_problem!` which stores the values of solution coefficients in vectors `η̂`, `η̇`, `ϕ̂`, `ϕ̇`, `ψ̂`, `ψ̇`. In our case only `η̂` will be further processed.

```@example 1
solve_problem!(η̂, η̇, ϕ̂, ϕ̇, ψ̂, ψ̇, β̂, β̇, p̂, κ, 𝒯, 𝒮, ℐ, M_s, M_b, Δt, O, N, χ, ξ, ζ, ℓ, d)
nothing # hide
```

Finally we are ready to plot evolution of absolute, real and imaginary values of a complex wave amplitude `η̂` over time `t` corresponding to one wave period.

```@example 1
t = range(start = 0, stop = N*Δt, step = Δt) # time vector
set_theme!(theme_latexfonts()) # use latex fonts
update_theme!(fontsize=10)
fig = Figure(size = (400, 300))
ax = Axis(fig[1, 1], xlabel = L"t/T", ylabel = L"4\hat{\eta}/H", xticks = 0:0.1:N, yticks = -1:0.5:1)
lines!(ax, t / T, 4 * abs.(η̂[1, O:end]) / H, label = L"|\hat{\eta}|")
lines!(ax, t / T, 4 * real.(η̂[1, O:end]) / H, label = L"Re(\hat{\eta})")
lines!(ax, t / T, 4 * imag.(η̂[1, O:end]) / H, label = L"Im(\hat{\eta})")
axislegend(ax, position = :lb)
fig
```

Now, if we want to plot a time-series of free-surface elevation at some location (e.g. corresponding to the middle of the fluid domain ``x_0 = \ell/2``) ``\eta(t[n], x_0)``, we use an `inverse_fourier_transform` function, like so

```@example 1
x_0 = ℓ/2
η(n) = inverse_fourier_transform(η̂[:, n], κ, x_0)
fig = Figure(size = (400, 300))
ax = Axis(fig[1, 1], xlabel = L"$t$ (s)", ylabel = L"$\eta$ (m)")
lines!(ax, t, η.(O:N+O), color = :blue, linewidth = 2)
limits!(ax, 0, T, -H, H)
fig
```

```@example 1
η(O), η(N+O), t
```