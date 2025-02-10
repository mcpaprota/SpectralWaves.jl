# Examples

## Regular linear waves

We are modelling linear and regular waves of length `L` and height `H` propagating in water of constant depth `d`. We consider one wave along the length of the domain `ℓ`. 

```@example 1
using SpectralWaves
using CairoMakie # plotting package

L = 2.0 # wavelength (m)
H = 0.1 # wave height (m)
d = 1.0 # water depth (m)
ℓ = L # fluid domain length (m) - one wave
nothing # hide
```

We need a wave period `T`.

```@example 1
k = 2π / L # wave number (rad/m)
ω = sqrt(g * k * tanh(k * d)) # angular wave frequency (rad/s)
T = 2π / ω # wave period (s)
nothing # hide
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

We define a function that calculates free-surface elevation at specified location `x` and time instant `n` and velocity components `u` and `w` at specified location (`x`, `z`) and time instant `n` using `water_surface`, `horizontal_velocity`, and `vertical_velocity` functions. Please note that the last argument to `water_velocity` specifies the axis of projection of velocity vector (1 - ``x``, 3 - ``z``)

```@example 1
η(x, n) = water_surface(p, x, n)
u(x, z, n) = water_velocity(x, z, n, 1)
w(x, z, n) = water_velocity(x, z, n, 3)
nothing # hide
```

Finally we are ready to plot results. Here, we plot the whole time series of free-surface elevation corresponding to the middle of a domain.

```@example 1
set_theme!(merge(theme_latexfonts(), Theme(fontsize=9))) # set latex fonts of size 9
fig = Figure(size = (400, 200)) 
ax = Axis(fig[1, 1], xlabel = L"t/T", ylabel = L"\eta/H")
lines!(ax, t / T, η.(ℓ/2, 1:length(t)) / H)
limits!(ax, 0, 1, -1, 1)
fig
```
